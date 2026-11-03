#!/usr/bin/env python3
"""Executa as células de código das aulas (.tex) e confere as saídas.

Uso:
    python3 scripts/verifica_codigos.py aulas/aula01.tex [--mat]
    python3 scripts/verifica_codigos.py --todas
    python3 scripts/verifica_codigos.py --notebooks     # gera notebooks/*.ipynb

Convenções no .tex (comentários LaTeX logo antes de \\begin{entrada}):
    % entradas: 3 | 4.5      -> respostas simuladas para input()
    % erro                   -> a célula deve gerar um erro (compara a última linha da saída)
    % nao-executar           -> a célula não é executada
"""
import ast
import contextlib
import io
import json
import re
import sys
import traceback
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
plt.show = lambda *a, **k: None   # no modo de verificação, show() não faz nada

RAIZ = Path(__file__).resolve().parent.parent
AMB = re.compile(r"\\begin\{(entrada|saida|saidamat)\}(?:\[[^\]]*\])?\n(.*?)\\end\{\1\}", re.S)
TITULO = re.compile(r"\\(chapter|section|subsection)\*?\{(.*)\}")


def blocos(texto):
    """Lista de (tipo, conteudo, marcas) na ordem em que aparecem."""
    out = []
    for m in AMB.finditer(texto):
        antes = texto[: m.start()].rstrip("\n").split("\n")[-3:]
        marcas = {}
        for linha in antes:
            linha = linha.strip()
            if linha.startswith("% entradas:"):
                marcas["entradas"] = [s.strip() for s in linha.split(":", 1)[1].split("|")]
            elif linha == "% erro":
                marcas["erro"] = True
            elif linha == "% nao-executar":
                marcas["pular"] = True
            elif linha.startswith("% figura:"):
                marcas["figura"] = linha.split(":", 1)[1].strip()
        out.append((m.group(1), m.group(2), marcas, m.start()))
    return out


def mostra(valor):
    try:
        import sympy
        if isinstance(valor, sympy.Basic) or (isinstance(valor, (list, tuple, dict, set)) and "sympy" in repr(type(valor))):
            return str(valor)
    except ImportError:
        pass
    return repr(valor)


def salva_figura(nome, valor):
    destino = RAIZ / "figuras" / f"{nome}.pdf"
    try:
        from sympy.plotting.plot import Plot
        if isinstance(valor, Plot):
            valor.save(str(destino))
            return
    except ImportError:
        pass
    if plt.get_fignums():
        plt.gcf().savefig(destino, bbox_inches="tight")
    else:
        print(f"  [aviso] célula sem figura para salvar: {nome}")


def executa(codigo, ns, entradas, figura=None):
    fila = list(entradas or [])

    def falso_input(msg=""):
        resp = fila.pop(0) if fila else "0"
        print(f"{msg}{resp}")
        return resp

    ns["input"] = falso_input
    arvore = ast.parse(codigo)
    ultimo = None
    if arvore.body and isinstance(arvore.body[-1], ast.Expr):
        ultimo = ast.Expression(arvore.body.pop().value)
    buf = io.StringIO()
    valor = None
    with contextlib.redirect_stdout(buf):
        exec(compile(arvore, "<celula>", "exec"), ns)
        if ultimo is not None:
            valor = eval(compile(ultimo, "<celula>", "eval"), ns)
            if valor is not None and not type(valor).__module__.startswith(("matplotlib", "sympy.plotting", "mpl_toolkits")):
                print(MARCA + mostra(valor))
    if figura:
        salva_figura(figura, valor)
    plt.close("all")
    return buf.getvalue(), valor


MARCA = "\x00VALOR\x00"


def normaliza(s):
    return "\n".join(l.rstrip() for l in s.strip("\n").split("\n")).rstrip()


def verifica(arquivo, mostrar_mat=False, preencher=False):
    texto = Path(arquivo).read_text(encoding="utf-8")
    trocas = []
    ns = {"__name__": "__main__"}
    lista = blocos(texto)
    problemas = 0
    ultimo_valor, ultima_saida, ultima_linha = None, "", 0
    for i, (tipo, conteudo, marcas, pos) in enumerate(lista):
        linha = texto[:pos].count("\n") + 1
        if tipo == "entrada":
            ultimo_valor, ultima_saida, ultima_linha = None, "", linha
            if marcas.get("pular"):
                continue
            try:
                ultima_saida, ultimo_valor = executa(conteudo, ns, marcas.get("entradas"), marcas.get("figura"))
                if marcas.get("erro"):
                    print(f"  [!] {arquivo}:{linha}: esperava erro e não houve")
                    problemas += 1
            except Exception as e:  # noqa: BLE001
                if marcas.get("erro"):
                    ultima_saida = f"{type(e).__name__}: {e}"
                else:
                    print(f"  [ERRO] {arquivo}:{linha}: {type(e).__name__}: {e}")
                    traceback.print_exc(limit=1)
                    problemas += 1
        elif tipo == "saida":
            if lista[i - 1][0] != "entrada" or lista[i - 1][2].get("pular"):
                continue
            esperado = normaliza(conteudo)
            saida_txt = ultima_saida
            if i + 1 < len(lista) and lista[i + 1][0] == "saidamat":
                saida_txt = saida_txt.split(MARCA)[0]       # o valor vai na caixa de fórmula
            obtido = normaliza(saida_txt.replace(MARCA, ""))
            if lista[i - 1][2].get("erro"):
                obtido = obtido.split("\n")[-1]
                esperado = esperado.split("\n")[-1]
            if esperado != obtido and preencher:
                m = AMB.search(texto, pos)
                trocas.append((m.start(2), m.end(2), obtido + "\n"))
                print(f"  [preenchido] {arquivo}:{ultima_linha}")
            elif esperado != obtido:
                problemas += 1
                print(f"  [DIFERENÇA] {arquivo}:{ultima_linha}")
                print("    --- no texto ---\n    " + esperado.replace("\n", "\n    "))
                print("    --- obtido ---\n    " + obtido.replace("\n", "\n    "))
        elif tipo == "saidamat" and preencher and "% manual" not in conteudo and ultimo_valor is not None:
            import sympy
            novo = "$\\displaystyle " + sympy.latex(ultimo_valor) + "$"
            if novo.strip() != conteudo.strip():
                m = AMB.search(texto, pos)
                trocas.append((m.start(2), m.end(2), novo + "\n"))
                print(f"  [mat preenchido] {arquivo}:{ultima_linha}")
        elif tipo == "saidamat" and mostrar_mat:
            try:
                import sympy
                lat = sympy.latex(ultimo_valor) if ultimo_valor is not None else ultima_saida
            except Exception:  # noqa: BLE001
                lat = repr(ultimo_valor)
            print(f"  [mat] {arquivo}:{ultima_linha}\n    texto : {conteudo.strip()}\n    sympy : {lat}")
    for ini, fim, novo in sorted(trocas, reverse=True):
        texto = texto[:ini] + novo + texto[fim:]
    if trocas:
        Path(arquivo).write_text(texto, encoding="utf-8")
    n = sum(1 for b in lista if b[0] == "entrada")
    print(f"{arquivo}: {n} células, {problemas} problema(s)")
    return problemas


# ---------------- geração de notebooks ----------------
def limpa_titulo(t):
    t = re.sub(r"\\texttt\{([^}]*)\}", r"`\1`", t)
    t = re.sub(r"\\py\{([^}]*)\}", r"`\1`", t)
    t = re.sub(r"\\(emph|textbf)\{([^}]*)\}", r"\2", t)
    t = t.replace("\\_", "_").replace("~", " ").replace("\\&", "&")
    return re.sub(r"\\[a-zA-Z]+\s*", "", t).replace("{", "").replace("}", "")


def notebook(arquivo, destino):
    texto = Path(arquivo).read_text(encoding="utf-8")
    eventos = []
    for m in TITULO.finditer(texto):
        nivel = {"chapter": "#", "section": "##", "subsection": "###"}[m.group(1)]
        eventos.append((m.start(), "markdown", f"{nivel} {limpa_titulo(m.group(2))}"))
    for m in re.finditer(r"\\begin\{exercicio\}(?:\[([^\]]*)\])?(?:\{[^}]*\})?(.*?)\\end\{exercicio\}", texto, re.S):
        corpo = m.group(2).replace("\\[", "$$").replace("\\]", "$$")
        partes = re.split(r"(\$\$.*?\$\$|\$[^$]*\$)", corpo, flags=re.S)
        for j in range(0, len(partes), 2):          # só o texto fora das fórmulas
            t = partes[j]
            t = re.sub(r"\\py\{([^}]*)\}|\\texttt\{([^}]*)\}", lambda k: "`" + (k.group(1) or k.group(2)) + "`", t)
            t = re.sub(r"\\item\s*", "\n- ", t)
            t = re.sub(r"\\(emph|textbf)\{([^}]*)\}", r"\2", t)
            t = t.replace("~", " ").replace("\\,", " ").replace("``", '"').replace("''", '"')
            t = re.sub(r"\\label\{[^}]*\}|\\ref\{[^}]*\}|\\(begin|end)\{[^}]*\}(\[[^\]]*\])?|\\[a-zA-Z]+\*?", "", t)
            partes[j] = t
        corpo = "".join(partes)
        eventos.append((m.start(), "markdown", "**Exercício.** " + " ".join(corpo.split()).replace(" - ", "\n- ")))
        eventos.append((m.start() + 1, "code", "# sua solução aqui\n"))
    for tipo, conteudo, marcas, pos in blocos(texto):
        if tipo == "entrada":
            eventos.append((pos, "code", conteudo.rstrip("\n")))
    eventos.sort(key=lambda e: e[0])
    celulas = []
    for _, tipo, fonte in eventos:
        cel = {"cell_type": tipo, "metadata": {}, "source": fonte.splitlines(keepends=True)}
        if tipo == "code":
            cel.update({"execution_count": None, "outputs": []})
        celulas.append(cel)
    nb = {"cells": celulas, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                          "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}
    Path(destino).write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")


if __name__ == "__main__":
    args = sys.argv[1:]
    aulas = sorted((RAIZ / "aulas").glob("aula*.tex"))
    if "--notebooks" in args:
        for a in aulas:
            notebook(a, RAIZ / "notebooks" / (a.stem + ".ipynb"))
            print("notebook:", a.stem + ".ipynb")
        sys.exit(0)
    alvo = aulas if "--todas" in args else [Path(a) for a in args if not a.startswith("--")]
    total = sum(verifica(a, "--mat" in args, "--preencher" in args) for a in alvo)
    sys.exit(1 if total else 0)
