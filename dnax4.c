//ate agora ja tem a comunicaçao entres os processos 
// p0 esta enviando toda a sequencia e o codon  para os outros processos
// o final da sequencia de um processo ja esta indo para o outro
// falta apenas comparar os cods com as seqs e depois guardar tudo em um arquivo



#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "arquivo.h"

int main(int argc, char** argv) {

	int my_rank, np;
	int tempo_logico = 0, tempo_logico_msg = 0;
	int sinal=1, pos = 0;// aqui diz se nao tem mais nenhum seq para enviar
	int posicao_cod=0;// onde o cod foi encontrado no seq
	int tamanho_seq[3];//tamanho da sequencia de DNA 
	int tamanho_cod[2];//tamanho do codon que queremos encontrar
	FILE *arq_seq = fopen("dna.in","r");// arquivo onde esta a sequencia de DNA
	FILE *arq_cod = fopen("query.in","r");
	FILE *arq_out = fopen("dna.out","a");
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);


	if(my_rank==0){
//______________________________________________________________________________________________________________________________________________
//______________________________________________________aqui estamos pegando o numero e o tamanho do maior cod__________________________________
//______________________________________________________e enviando para os outros processos_____________________________________________________
		fseek(arq_cod,0,SEEK_SET);
		tamanho_cod[0]=ncod(arq_cod);
		fseek(arq_cod,0,SEEK_SET);
		tamanho_cod[1]=maior_cod(arq_cod);
			
		
		for (int i=1;i<np;i++){
			// enviar aqui
			tempo_logico++;
			MPI_Send(&tempo_logico, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
			MPI_Send(tamanho_cod,2,MPI_INT,i,0,MPI_COMM_WORLD);
		
		}

//______________________________________________________________________________________________________________________________________________
//________________________________________agora vamos dividir a sequencia e enviara a parte de cada um__________________________________________
		int resto_divisao=0;
		
		int verifica=0;
		tamanho_seq[0]=0;
		tamanho_seq[1]=1;
		fseek(arq_seq,0,SEEK_SET);
		
		while(sinal==1){
			verifica=TamanhoSequencia(arq_seq);
			if(verifica==-1){
				sinal=0;
			}
			// esse if é para poder enviar as outras sequencias 	
			if(sinal!=0){
				tamanho_seq[0]=tamanho_seq[0]+verifica;
				tamanho_seq[1]=TamanhoSequencia(arq_seq)-1;
				//printf("tamanho_seq[0]:%d \n",tamanho_seq[0]);
				resto_divisao=tamanho_seq[1]%(np-1);
				
				//tamanho_seq[2]=resto_divisao;//passa o resto da divisao
				
				tamanho_seq[1]=tamanho_seq[1]/(np-1);
				//printf("tamanho_seq[1]:%d \n",tamanho_seq[1]);
				//printf("resto_divisao:%d \n",resto_divisao);
				//printf("_____________________________________________________\n");
				char c[tamanho_seq[1]];
				// estamos enviando a parte que cada um tem que pegar no arquivo
				tamanho_seq[2]=0;
				for(int i =1;i<np;i++){
				
					if(i<=resto_divisao){			
						tamanho_seq[1]=tamanho_seq[1]+1;
						//enviar aqui
						tempo_logico++;
						MPI_Send(&tempo_logico, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
						MPI_Send(tamanho_seq,3,MPI_INT,i,0,MPI_COMM_WORLD);
						tamanho_seq[2]=tamanho_seq[2]+tamanho_seq[1];
						tamanho_seq[1]=tamanho_seq[1]-1;
						tamanho_seq[0]=tamanho_seq[0]+tamanho_seq[1]+1;
					}
					else{
						//enviar aqui
						tempo_logico++;
						MPI_Send(&tempo_logico, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
						MPI_Send(tamanho_seq,3,MPI_INT,i,0,MPI_COMM_WORLD);
						tamanho_seq[2]=tamanho_seq[2]+tamanho_seq[1];
						tamanho_seq[0]=tamanho_seq[0]+tamanho_seq[1];			
					}
				}
			}
			tamanho_seq[0]=tamanho_seq[0]+1;
			//enviar aqui
			tempo_logico++;
			MPI_Send(&tempo_logico, 1, MPI_INT, 1, 3, MPI_COMM_WORLD);
			MPI_Send(&sinal,1,MPI_INT,1,1,MPI_COMM_WORLD);
		}
		
		
	//__________________________________________________________________________________________________________________________________________
	}
	else{//my_rank!=0	
	//________________________________________________________________________________________________________________________________________
	//___________________________________________estamos recebendo o tamanho e o numero de cod _________________________________________________
	//___________________________________________e criando uma matriz com os cods_______________________________________________________________
		char **Matriz_cod;
		//receber aqui
		MPI_Recv(&tempo_logico_msg, 1, MPI_INT, 0,3, MPI_COMM_WORLD, &status);
		tempo_logico = max(tempo_logico, tempo_logico_msg);
		printf("Eu, processo %d, recebi %d no tempo lógico.\n", my_rank, tempo_logico);
		MPI_Recv(tamanho_cod,2,MPI_INT,0,0,MPI_COMM_WORLD,&status);
		Matriz_cod=alocarMatriz(tamanho_cod[0],tamanho_cod[1],arq_cod);
		
	//_____________________________________________________________________________________________________________________________________________	
	//_____________________________________________aqui recebemos os parametro do seq _____________________________________________________
	//_____________________________________________geramos o seq e enviamos a parte final para o proximo processo_________________________		
		if(my_rank==1){
			// esse while é para poder receber as outras seq
			int genoma = 1;
			while(sinal==1){
			// recebendo parametros do rank=0
				//receber aqui
				MPI_Recv(&tempo_logico_msg,1,MPI_INT,0,3,MPI_COMM_WORLD,&status);
				tempo_logico = max(tempo_logico, tempo_logico_msg);
				MPI_Recv(&sinal,1,MPI_INT,0,1,MPI_COMM_WORLD,&status);
				printf("Eu, processo %d, recebi %d no tempo lógico.\n", my_rank, tempo_logico);

				tempo_logico++;
				MPI_Send(&tempo_logico, 1, MPI_INT,my_rank+1,3,MPI_COMM_WORLD);
				MPI_Send(&sinal,1,MPI_INT,my_rank+1,1,MPI_COMM_WORLD);

				if(sinal==1){
					//receber aqui
					MPI_Recv(&tempo_logico_msg,1,MPI_INT,0,3,MPI_COMM_WORLD,&status);
					tempo_logico = max(tempo_logico, tempo_logico_msg);
					MPI_Recv(tamanho_seq,3,MPI_INT,0,0,MPI_COMM_WORLD,&status);
					printf("Eu, processo %d, recebi %d no tempo lógico.\n", my_rank, tempo_logico);
					char c[tamanho_seq[1]+1];
					
					// gera o seq
					pega_seq(tamanho_seq[0], tamanho_seq[1],arq_seq,c);
					
					c[tamanho_seq[1]]='\0';
					
					
					// envia parte final do seq
					//receber aqui
					tempo_logico++;
					MPI_Send(&tempo_logico, 1, MPI_INT,2,3,MPI_COMM_WORLD);
					MPI_Send(&c[tamanho_seq[1]-tamanho_cod[1]],tamanho_cod[1],MPI_BYTE,2,0,MPI_COMM_WORLD);

					//enviar aqui
					tempo_logico_msg = saida(Matriz_cod, c, tamanho_cod[0], my_rank, arq_out, genoma, tamanho_seq[1]*(my_rank-1), np, status,
										tempo_logico);
					tempo_logico = max(tempo_logico, tempo_logico_msg);
					genoma++;
					
					
				}
				
				//printf("sinal:%d\n",sinal);
				//printf("string ------------------- %s\n", x);
			}
	//_______________________________________________________________________________________________________________________________________
		}
		else{
	//____________________________________________________________________________________________________________________________________
	// ___________________________________________recebendo final do seq do processo anterios______________________________________________
	//_____________________________________________aqui recebemos os parametro do seq _____________________________________________________
	//_____________________________________________geramos o seq e enviamos a parte final para o proximo processo_________________________
			// esse while é para poder receber as outras seq
			int genoma = 1;
			
			while(sinal==1){
				//receber aqui
				MPI_Recv(&tempo_logico_msg,1,MPI_INT,my_rank-1,3,MPI_COMM_WORLD,&status);
				tempo_logico = max(tempo_logico, tempo_logico_msg);
				MPI_Recv(&sinal,1,MPI_INT,my_rank-1,1,MPI_COMM_WORLD,&status);
				printf("Eu, processo %d, recebi %d no tempo lógico.\n", my_rank, tempo_logico);

				if(my_rank!=(np-1)){
					//enviar aqui
					tempo_logico++;
					MPI_Send(&tempo_logico, 1, MPI_INT,my_rank+1,3,MPI_COMM_WORLD);
					MPI_Send(&sinal,1,MPI_INT,my_rank+1,1,MPI_COMM_WORLD);
				}
				if(sinal==1){
					// recebendo parametros do rank=0
					//receber aqui
					MPI_Recv(&tempo_logico_msg,1,MPI_INT,0,3,MPI_COMM_WORLD,&status);
					tempo_logico = max(tempo_logico, tempo_logico_msg);
					MPI_Recv(tamanho_seq,3,MPI_INT,0,0,MPI_COMM_WORLD,&status);
					printf("Eu, processo %d, recebi %d no tempo lógico.\n", my_rank, tempo_logico);
					
					
					char c[tamanho_seq[1]+1+tamanho_cod[1]];
					
					// recebendo final do seq do processo anterios
					//receber aqui
					MPI_Recv(&tempo_logico_msg,1,MPI_INT,my_rank-1,3,MPI_COMM_WORLD,&status);
					tempo_logico = max(tempo_logico, tempo_logico_msg);
					MPI_Recv(c,tamanho_cod[1],MPI_INT,my_rank-1,0,MPI_COMM_WORLD,&status);
					printf("Eu, processo %d, recebi %d no tempo lógico.\n", my_rank, tempo_logico);
					
					// gera o seq
					pega_seq(tamanho_seq[0], tamanho_seq[1],arq_seq,&c[tamanho_cod[1]] );
					// envia parte final do seq obs o ultimo nao envia o seq
					if(my_rank!=(np-1)){
						//enviar aqui
						tempo_logico++;
						MPI_Send(&tempo_logico, 1, MPI_INT,my_rank+1,3,MPI_COMM_WORLD);
						MPI_Send(&c[tamanho_seq[1]],tamanho_cod[1],MPI_BYTE,my_rank+1,0,MPI_COMM_WORLD);
					}
					//printf("seeq - > %s \n",c);
					int mult = 1;
					mult = tamanho_seq[2]-tamanho_cod[1];
					//enviar aqui
					tempo_logico_msg = saida(Matriz_cod, c, tamanho_cod[0], my_rank, arq_out, genoma, mult, np, status, tempo_logico);
					tempo_logico = max(tempo_logico, tempo_logico_msg);
					genoma++;
						
				}
			}
	//_______________________________________________________________________________________________________________________________________
		}	


	fflush(stdout);
	}
	
	fclose(arq_seq);
	fclose(arq_cod);
	fclose(arq_out);
	MPI_Finalize();
}

























