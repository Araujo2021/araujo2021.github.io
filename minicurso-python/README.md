# Python para Matemática — material de 10 aulas (MAT 172)

Material em LaTeX para a parte final da disciplina MAT 172, com Python no JupyterLab.
Baseado no livro *Python para Matemáticos* (A. L. L. Souza, SBM, 2023, cap. 3 a 6)
e no *Tutorial rápido de Python para Matemáticos* (R. M. Martins, 2022).

## Arquivos prontos

| Arquivo | Conteúdo |
|---|---|
| `curso-professor.pdf` | versão do professor: texto completo + notas ao professor + soluções de todos os exercícios |
| `curso-aluno.pdf` | versão do estudante: mesmo texto, sem notas e sem soluções |
| `notebooks/aulaNN.ipynb` | um notebook por aula com os títulos, todas as células de código e os enunciados dos exercícios, para abrir no JupyterLab |

## Plano das aulas

| Aula | Tema |
|---|---|
| 1 | Revisão de Python I: `print`, `input`, tipos, operações, lógica, `if`/`elif`/`else`, `for`, `while`, indentação |
| 2 | Revisão de Python II: listas, tuplas, dicionários, funções (`def`, `return`, `lambda`) |
| 3 | Módulo `math`; raiz quadrada pelo método babilônico; primalidade; bissecção |
| 4 | NumPy I: arrays, operações vetoriais, máscaras, somas, regra dos trapézios, Monte Carlo |
| 5 | NumPy II: matrizes, `linalg`, eliminação de Gauss, autovalores, polinômios e ajuste |
| 6 | Gráficos com Matplotlib (2D, paramétricas, polares, curvas de nível, 3D) |
| 7 | SymPy I: expressões, simplificação, equações, sistemas, somatórios |
| 8 | SymPy II: limites, derivadas, tangentes, máximos e mínimos, Taylor |
| 9 | SymPy III: integrais, Riemann, áreas, volumes, integrais duplas |
| 10 | SymPy IV: matrizes simbólicas, EDOs (`dsolve`), campo de direções; atividade final |

Apêndices: tabela Maxima → Python, referência rápida, soluções (só na versão do professor).

## Estrutura

```
principal.tex          corpo do documento (inclui aulas e apêndices)
curso-professor.tex    gera a versão do professor
curso-aluno.tex        gera a versão do estudante
preambulo.tex          pacotes, cores e caixas (células no estilo Jupyter, Atenção, Dica, Maxima...)
apresentacao.tex       apresentação e plano geral
aulas/aula01..10.tex   uma aula por arquivo
apendices/             Maxima→Python, referência rápida, soluções, referências
figuras/               figuras (geradas pelos scripts)
scripts/               verificação dos códigos, geração de figuras e notebooks
```

## Como compilar

```bash
latexmk -pdf curso-professor.tex
latexmk -pdf curso-aluno.tex
```

## Como os códigos foram conferidos

Todas as células do texto são executadas por `scripts/verifica_codigos.py`, que compara a saída
mostrada no PDF com a saída real do Python (Python 3.13, NumPy 2.5, SymPy 1.14, Matplotlib 3.11).

```bash
python3 scripts/verifica_codigos.py --todas                 # confere todas as aulas
python3 scripts/verifica_codigos.py apendices/solucoes.tex  # confere as soluções
python3 scripts/verifica_codigos.py aulas/aula07.tex --preencher   # atualiza as saídas no .tex
python3 scripts/gera_figuras.py                             # figuras ilustrativas
python3 scripts/verifica_codigos.py --notebooks             # gera notebooks/*.ipynb
```

Convenções nos `.tex`: comentários `% entradas: 3 | 4` (respostas simuladas para `input()`),
`% erro` (célula que deve gerar erro) e `% figura: nome` (salva o gráfico da célula em `figuras/nome.pdf`).
Ao editar um código, rode o verificador com `--preencher` para atualizar a saída mostrada.
