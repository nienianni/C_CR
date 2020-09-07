#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
int NETWORK_SIZE_ini = 30;//用于搭建网络
double PROBABILITY_OF_EAGE = 0.8;
int NETWORK_SIZE = 50;
int e_link_num=20;//每次新加入节点连接的边数
int ** adjacentMatrix_w;
int *sum_degree;//记录每个节点的度数
double ** adjacentMatrix_R;

double pd;
double lambda;
double mu = 0.2;
double la ;
int tmax = 401;
int t_test = 1;


FILE * fp1;
FILE * fp2;
FILE *fp3;
FILE *fp4;
int *city_people_num;
double *sum_w;
double ** peopleTodes_num;
double*I_pro_ini;
double *P;
double *PI;
double *rho;
double *rho_t;
double ** adjacentMatrix_w_alpha;
double** contact_M;
double *sum;
struct people
{
	int sta;
	char ini_state;
	int des;
	char fin_state;
}** city;

void load_city_people_num();
void initial();
//void generateNetwork_ini();
//void generate_SF_Network_ini();//生成SF网络
void load_Matrix_w();//将邻接矩阵变为带权的
void load_Matrix_R_fenzi();
void load_Matrix_R();
void load_struct_people();
void load_peopleTodes_num();
void load_I_pro_ini();
void load_P();
void load_PI();
void load_rho(int t);
void load_rho_t(int t);
void load_limit_rho();
double load_I_pro_zong();
void load_write_file();
void load_Matrix_w_alpha();
void load_contact_M();
void load_write_M();

int main()
{
	int i, j;
	load_city_people_num();
	initial();
	load_struct_people();//y
	load_I_pro_ini();//y
	load_Matrix_w();
	
	load_write_file();
	//load_write_M();
	return 0;
}
void load_city_people_num()
{
	city_people_num = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	int i, a;
	errno_t err;
	err = fopen_s(&fp1, "peoplenum_a.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		fscanf_s(fp1, "%d ", &a);
		city_people_num[i] = a;
	}
	fclose(fp1);
}
void initial()
{
	adjacentMatrix_w = (int**)malloc(sizeof(int *) * NETWORK_SIZE); //分配指针数组
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_w[i] = (int *)malloc(sizeof(int) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	adjacentMatrix_R = (double**)malloc(sizeof(double *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_R[i] = (double*)malloc(sizeof(double) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	city = (struct people**)malloc(sizeof(struct people *) * NETWORK_SIZE);//有多少个城市,city为指向指针的指针
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		city[i] = (struct people*)malloc(sizeof(struct people)* city_people_num[i]);//叶子城市里有多少个人
	}
	sum_w = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	peopleTodes_num = (double**)malloc(sizeof(double*) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		peopleTodes_num[i] = (double*)malloc(sizeof(double) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	I_pro_ini = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	P = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	PI = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	rho = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	rho_t = (double *)malloc(sizeof(double)*t_test);
	contact_M = (double**)malloc(sizeof(double *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		contact_M[i] = (double *)malloc(sizeof(double) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	adjacentMatrix_w_alpha = (double**)malloc(sizeof(double *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_w_alpha[i] = (double *)malloc(sizeof(double) * NETWORK_SIZE);//分配每个指针指向的数组
	}

	sum = (double *)malloc(sizeof(double)*NETWORK_SIZE);
}
void load_struct_people()
{
	int i, j, m;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < city_people_num[i]; j++)
		{
			city[i][j].sta = i;
			if (i == 0 && j < 10)
			{
				city[i][j].ini_state = 'I';
			}
			else
			{
				city[i][j].ini_state = 'S';
			}
			city[i][j].des = 100;
			city[i][j].fin_state = 'N';
		}
	}
}
void load_I_pro_ini()
{
	int i, j, m, S, I;
	for (m = 0; m < NETWORK_SIZE; m++)
	{
		S = 0;//记录每个种群的S人数
		I = 0;//记录每个种群的I人数
		for (i = 0; i < NETWORK_SIZE; i++)
		{
			for (j = 0; j < city_people_num[i]; j++)
			{
				if (city[i][j].sta == m && city[i][j].ini_state == 'S')
					S++;
				if (city[i][j].sta == m && city[i][j].ini_state == 'I')
					I++;
			}
		}
		I_pro_ini[m] = (double)I / (I + S);
	}
}
void load_Matrix_w()
{
	int i, j;
	int a;
	errno_t err;
	err = fopen_s(&fp2, "w_a.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			fscanf_s(fp2, "%d ", &a);
			adjacentMatrix_w[i][j] = a;
		}
	}
	fclose(fp2);
	/*for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] != 0)
			{
				adjacentMatrix_w[i][j] = 10;
			}
		}
	}*/
}

void load_Matrix_R_fenzi()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] == 0)
			{
				adjacentMatrix_R[i][j] = 0;
			}
			else
			{
				adjacentMatrix_R[i][j] = pow(adjacentMatrix_w[i][j] * (rho[j]+0.0000000001), la);// 
			}
			/*if (adjacentMatrix_w[i][j] == 0 || rho[j] < 0.000000000000000000001)
			{
				adjacentMatrix_R[i][j] = 0.00001;
			}
			else
			{
				adjacentMatrix_R[i][j] = pow(adjacentMatrix_w[i][j] * rho[j], la)+0.00001;
			}*/
		}
	}
}
void load_Matrix_R()
{
	int i, j;
	double s;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		s = 0;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			s = adjacentMatrix_R[i][j] + s;
		}
		sum_w[i] = s;
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			adjacentMatrix_R[i][j] = adjacentMatrix_R[i][j] / sum_w[i];
		}
	}
}
void load_peopleTodes_num()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (i == j)
			{
				peopleTodes_num[j][i] = (1 - pd)*city_people_num[i];
			}
			if (i != j)
			{
				peopleTodes_num[j][i] = pd * adjacentMatrix_R[j][i] * city_people_num[j];
			}
		}
	}
}
void load_P()
{
	int i, j;
	double a;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		a = 1;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			a = (pow((1 - lambda * rho[j]), peopleTodes_num[j][i]))*a;
		}
		P[i] = 1 - a;
	}
}
void load_PI()
{
	int i, j;
	double sum;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		sum = 0;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			sum = (adjacentMatrix_R[i][j] * P[j]) + sum;
		}
		PI[i] = (1 - pd)*P[i] + pd * sum;
	}
}
void load_rho(int t)
{
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		if (t == 0)
			rho[i] = I_pro_ini[i];
		else
		{
			rho[i] = (1 - mu)*rho[i] + (1 - rho[i])*PI[i];
			/*if (rho[i] < 0.0000000000000001&&rho[i]!=0)
			{
				printf("%.30If ", rho[i]);
			}*/
		}
	}
}
void load_rho_t(int t)
{
	int i;
	double sum = 0;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		sum = sum + rho[i];
	}
	rho_t[t - (tmax - t_test)] = (double)sum / NETWORK_SIZE;
}
void load_limit_rho()
{
	int t, i, j;
	for (t = 0; t < tmax; t++)
	{
		//printf("la:%If pd:%If lambda:%If t:%d \n", la, pd,lambda,t);
		load_rho(t);
		/*printf("rho\n");
		for (i = 0; i < NETWORK_SIZE; i++)
		{
			printf("%lf ", rho[i]);

		}
		printf("\n");
		//printf("\n");*/
		if (t >= (tmax - t_test))
		{
			load_rho_t(t);
		}
		/*printf("W:\n");
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%d ", adjacentMatrix_w[j][i]);
			}
			printf("\n");
		}
		printf("\n");*/
		load_Matrix_R_fenzi();
		/*printf("R_fenzi:\n");
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%f ", adjacentMatrix_R[j][i]);
				if (adjacentMatrix_R[j][i] < 0)
				{
					printf("%f ", adjacentMatrix_R[j][i]);
					printf("%d,%d\n", j, i);
				}

			}
			printf("\n");
		}
		printf("\n");*/
		load_Matrix_R();
		/*printf("R:\n");
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%f ", adjacentMatrix_R[j][i]);
			}
		}
		printf("\n");*/
		load_peopleTodes_num();
		/*printf("peopletodes:\n");
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%f ", peopleTodes_num[j][i]);
			}
			printf("\n");
		}
		printf("\n");*/
		load_P();
		load_PI();
	}
}

double load_I_pro_zong()
{
	int i;
	double pro_zong, zong = 0;
	/*for (i = 0; i < t_test; i++)
	{
		printf("%d:%f\n",i, rho_t[i]);
	}*/
	for (i = 0; i < t_test; i++)
	{
		zong = rho_t[i] + zong;
	}
	pro_zong = (double)zong / t_test;
	return pro_zong;
}
void load_write_file()
{
	double pro;
	int i, j;
	errno_t err;
	err = fopen_s(&fp3, "I_pd_L_a_3_2.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (la = -1; la <= 1; la = la + 1)//11
	{
		fprintf_s(fp3, "la=%f\n", la);
		for(lambda = 0.004; lambda <= 0.009; lambda = lambda + 0.002)//50
		{
			for (pd = 0; pd <= 1.01; pd = pd + 0.02)//50
			{
				load_limit_rho();
				pro = load_I_pro_zong();
				fprintf_s(fp3, "%f ", pro);
				//printf("%f %f %f %.2f", pd, la, lambda, pro);
			}
			fprintf_s(fp3, "\n");
			//printf("\n");
		}
		fprintf_s(fp3, "\n");
		fprintf_s(fp3, "\n");
		fprintf_s(fp3, "\n");
	}
	fclose(fp3);
}
/*void load_write_file()
{
	double pro;
	int i, j;
	errno_t err;
	err = fopen_s(&fp3, "help1.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (pd = 0.3; pd <= 0.8; pd = pd + 0.2)//11
	{
		fprintf_s(fp3, "pd=%f\n", pd);
		for (la = -6; la <= 6; la = la + 3)//50
		{
			for (lambda = 0; lambda <= 0.01; lambda = lambda + 0.0002)//50
			{
				load_limit_rho();
				pro = load_I_pro_zong();
				fprintf_s(fp3, "%f ", pro);
				//printf("%f %f %f %.2f", pd, la, lambda, pro);
			}
			fprintf_s(fp3, "\n");
			//printf("\n");
		}
		fprintf_s(fp3, "\n");
		fprintf_s(fp3, "\n");
		fprintf_s(fp3, "\n");
	}
	fclose(fp3);
}*/

void load_Matrix_w_alpha()
{
	int i, j;
	double s;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		s = 0;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] != 0)
			{
				s = s + pow(adjacentMatrix_w[i][j], la);
			}
		}
		sum[i] = s;//存储每一行的和，共i行
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] == 0)
			{
				adjacentMatrix_w_alpha[i][j] = 0;
			}
			else
				adjacentMatrix_w_alpha[i][j] = pow(adjacentMatrix_w[i][j], la) / sum[i];//出现误差
		}
	}
}
void load_contact_M()
{
	double sum;
	int i, j, l, m, n;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			sum = 0;
			for (l = 0; l < NETWORK_SIZE; l++)
			{
				sum = sum + (adjacentMatrix_w_alpha[i][l] * adjacentMatrix_w_alpha[j][l]);
			}
			if (i == j)
			{
				contact_M[i][j] = (1 - pd)*(1 - pd) *city_people_num[j] + pd * (1 - pd)*city_people_num[j] * (adjacentMatrix_w_alpha[i][j] + adjacentMatrix_w_alpha[j][i]) + pd * pd*city_people_num[j] * sum;
			}
			if (i != j)
			{
				contact_M[i][j] = pd * (1 - pd)*city_people_num[j] * (adjacentMatrix_w_alpha[i][j] + adjacentMatrix_w_alpha[j][i]) + pd * pd*city_people_num[j] * sum;
			}
		}
	}
}
void load_write_M()
{
	double pro;
	int i, j;
	errno_t err;
	err = fopen_s(&fp4, "M_a_3.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (la = -1; la <= 1; la = la + 1)
	{
		fprintf_s(fp4, "la=%f\n", la);
		for (pd = 0; pd <= 1.01; pd = pd + 0.02)
		{
			load_Matrix_w_alpha();
			load_contact_M();
			//fprintf_s(fp3, "pd=%f\n", pd);
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				for (j = 0; j < NETWORK_SIZE; j++)
				{
					fprintf_s(fp4, "%f ", contact_M[i][j]);
				}
				fprintf_s(fp4, "\n");
			}
		}
		fprintf_s(fp4, "\n");
		fprintf_s(fp4, "\n");
		fprintf_s(fp4, "\n");
	}
	fclose(fp4);
}