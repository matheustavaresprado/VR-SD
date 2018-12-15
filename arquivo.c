#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int TamanhoSequencia(FILE *arq){
	char c;
	
	int contador=0;
	c=fgetc(arq);
	if(c==EOF){
		return -1;
	}
	while(c!=EOF && c!='\n'){

		contador++;
		c=fgetc(arq);

	}
	contador++;
	return contador;
	
}


int ncod(FILE *arq){
	fseek(arq,0,SEEK_SET);
	int contador=0;
	while(TamanhoSequencia(arq)!=-1){
		contador++;
	}
	fseek(arq,0,SEEK_SET);
	return (contador/2);	
}

int maior_cod(FILE *arq){
	fseek(arq,0,SEEK_SET);
	int maior=0;
	TamanhoSequencia(arq);
	int atual=TamanhoSequencia(arq);
	while(atual!=-1){
		if(atual>maior){
			maior = atual;		
		}
		TamanhoSequencia(arq);
		atual=TamanhoSequencia(arq);
	}
	fseek(arq,0,SEEK_SET);
	return maior-1;

}




char** alocarMatriz(int Linhas,int Colunas,FILE *arq){ //Recebe a quantidade de Linhas e Colunas como Parâmetro
 
	  int i,j; //Variáveis Auxiliares
	  Colunas= Colunas+1;
	  char **m = (char**)malloc(Linhas * sizeof(char*)); //Aloca um Vetor de Ponteiros
	 
	  for (i = 0; i < Linhas; i++){ //Percorre as linhas do Vetor de Ponteiros
		   m[i] = (char*) malloc(Colunas * sizeof(char)); //Aloca um Vetor de Inteiros para cada posição do Vetor de Ponteiros.
		   for (j = 0; j < Colunas; j++){ //Percorre o Vetor de Inteiros atual.
				m[i][j] = '\0'; //Inicializa com \0.
		   }
	 }

	char c;
	for(i=0;i<Linhas;i++){
		TamanhoSequencia(arq);
		for(j=0;j<Colunas;j++){
			c=fgetc(arq);
			if(c=='\n' || c==EOF){//o sistema estava colocando caracteres estranhos no ultimo campo, botei o eof tbm;
				m[i][j]='\0';
				j=Colunas;
			}
			else{
				m[i][j]=c;
			}

		}
	}
	return m; //Retorna o Ponteiro para a Matriz Alocada
}

int Compara_STR(char *cadeia,char *codon){
	int p=0;
	for(int i=0;cadeia[i]!='\0';i++){
		if(cadeia[i]==codon[p]){
			p++;
		}
		else{	
			p=0;
		}
		if(codon[p]==10 || codon[p]==13 || codon[p] == '\0' || codon[p] == '\n' || codon[p] == EOF){
			codon[p]='\0';
			return i-(p-2);
		}
		
	}
	fflush(stdout);
	return -1;
}

void pega_seq(int inicio, int fim, FILE *arq, char *a){

	fseek(arq,0,SEEK_SET);
	fseek(arq,inicio,SEEK_SET);
	int i;
	for( i = 0;(i<fim || *(a+i)=='\n'); i++){
		*(a+i)=fgetc(arq);
	}
	*(a+i)='\0';
}

int max(int n1, int n2){
	if(n1>n2)return n1;
	return n2;
}


int saida(char **querys, char *seq, int n, int nu_proc, FILE *out, int genoma, int pos_init, int np, MPI_Status status, int tempo_logico){
	int result = -1;
	int sinal = -1;
	int tempo_logico_msg = tempo_logico;
	for (int i = 0; i<n; i++){
		
		result = Compara_STR(seq, &querys[i][0]);
		if(nu_proc == 1){
			//enviar aqui
			tempo_logico++;
			MPI_Send(&tempo_logico, 1, MPI_INT,nu_proc+1,3,MPI_COMM_WORLD);
			MPI_Send(&result,1,MPI_INT,nu_proc+1,8,MPI_COMM_WORLD);
		}
		if(nu_proc != 1){
			//receber aqui
			MPI_Recv(&tempo_logico_msg,1,MPI_INT,nu_proc-1,3,MPI_COMM_WORLD,&status);
			printf("Processo: %d | Remetente: %d | Tempo: %d | Tempo recebido: %d\n", nu_proc, nu_proc-1, tempo_logico, tempo_logico_msg);
			tempo_logico = max(tempo_logico, tempo_logico_msg);
			MPI_Recv(&sinal,1,MPI_INT,nu_proc-1,8,MPI_COMM_WORLD,&status);
			
			if(sinal != -1){
				result = 1;
			}
		}
		if(nu_proc != (np-1) && nu_proc !=1){
			//enviar aqui
			tempo_logico++;
			MPI_Send(&tempo_logico, 1, MPI_INT,nu_proc+1,3,MPI_COMM_WORLD);
			MPI_Send(&result,1,MPI_INT,nu_proc+1,8,MPI_COMM_WORLD);
		}
		if(result!=-1 && (-1==sinal)){//achou		
			fseek(out,0,SEEK_END);
			fprintf(out,"Query string %d\n>Escherichia coli partial genome (%d)\n%d ->%s %d \n", (i+1), genoma, (result+pos_init-1), &querys[i][0],nu_proc);
		}
		else{
			//fprintf(out,"Query string %d\nNOT FOUND\n", i);
		}
	}
	return tempo_logico;
}






















