# Purkinje-Microscopic-MEF

Resolve a equação do monodomínio utilizando o Método dos Elementos Finitos usando Elementos de Hermite. O projeto está dividido em 3 partes.

  - Skeleton_Mesh
  - Mesh_Generator
  - Solver

# Skeleton_Mesh

  - Constrói o esqueleto da malha a ser resolvida.
  - Atualmente está voltada para as simulações ligadas ao estudo do problema de Fonte-Sumidouro nas bifurcações das fibras.
  - A malha de saída do programa está configurada com a extensão .vtk para poder ser visualizada no Paraview. 
  - Para usar basta passar como argumento para o programa o número de bifurcações e o tamanho de cada fibra.

```sh
$ make
$ ./skeletonMesh <xMax> <num_bif> <out_VTK_file>
```

# Mesh_Generator

  - Constrói a malha de elementos finitos a partir da malha vinda do Skeleton_Mesh.
  - Para usar basta passar como argumento para o programa o número de elementos que cada fibra irá possuir e o tamanho de cada fibra (deve ser o mesmo que o usado no Skeleton_Mesh).

```sh
$ make
$ ./purkinjeMiocardiumMEF <xMax> <nElem> <VTK_file>
```

# Solver

  - Resolve a equação monodomínio utilizando o MEF para geração do potencial de ação dos miócitos da fibra.
  - Considera conservação de corrente de acordo com a Lei de Kirchoff nas bifurcações.
  - A EDP associada do problema é resolvida usando decomposição LU com pivoteamento. A cada passo de tempo se faz retro+pos substituições.
  - A sistema não linear de EDOs associado é resolvido usando Euler Explícito.
  - Atualmente está versão está configurada para o modelo celular de Noble.
  - Para usar basta passar como argumento para o programa o passo de tempo 'dt', o período máximo da simulação 't_max' e o arquivo da malha gerado a partir do Mesh_Generator. 
  - Solução fica armazena na pasta VTK contendo os valores do potencial transmembrânico e da corrente de todos os miócitos.
  - Para visualizar a simulação abrir os arquivos da pasta VTK no Paraview.

```sh
$ make
$ ./purkinjeFEM <dt> <t_max> <mesh_file>
```
