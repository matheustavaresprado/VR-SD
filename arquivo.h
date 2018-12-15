#include <stdio.h>
#include <stdlib.h>
#include "arquivo.c"
#include <mpi.h>

// essa funçao foi feita para contar quantos cod estao no arquivo "query.txt" 
// ela apenas recebe um arquivo e retorna a metade da quantidade de linhas  
int ncod(FILE *arq);
//essa funçao recebe um arquivo e devolve um int com o numero de caracteres ate o primeiro \n ou seja, ate pular lina
// caso o arquivo esteja no final ela devolve -1
int TamanhoSequencia(FILE *arq);
//ela recebe dois ponteiros para char e procura o char codon dentro do Cadeia caso ache ela retora a posiçao que ele foi encontrado
int Compara_STR(char *cadeia,char *codon);
// aqui ele olha os cod e devolve o numero de caracteres do maior
int maior_cod(FILE *arq);
// essa é usada para criar uma matriz de cod
char** alocarMatriz(int Linhas,int Colunas,FILE *arq);



int max(int n1, int n2);
void pega_seq(int inicio, int fim,FILE *arq,char *a );
int saida(char **query, char *seq, int n, int nu_proc,FILE *out, int genoma, int pos_init, int np, MPI_Status status, int tempo_logico);

