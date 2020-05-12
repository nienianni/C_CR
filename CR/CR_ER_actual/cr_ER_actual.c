#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

int NETWORK_SIZE = 50;
double PROBABILITY_OF_EAGE = 0.5;
float ** adjacentMatrix_w;
double ** adjacentMatrix_R;
double pd;
double lambda;
double mu = 0.2;
double la=-2;

//FILE * fp3;
FILE * fp4;
FILE *fp;
FILE * fp5;
int *city_people_num;
double *sum_w;
int ** peopleTodes_num;
double *I_pro_ini;
double *P;
double *PI;
double *rho;
struct people
{
	int sta;
	char ini_state;
	int des;
	char fin_state;
}** city;

void load_city_people_num();
void initial();
void generateNetwork_ini();
void load_Matrix_w();//将邻接矩阵变为带权的
void load_Matrix_R_fenzi();
void load_Matrix_R();
void load_struct_people();
void load_peopleTodes_num();
void load_I_pro_ini();
void load_P();
void load_PI();
void load_rho(int t);//100个步长后的每个集群的感染占比
void load_limit_rho();
double load_I_pro_zong();
void load_write_file();

int main()
{
	int i, j;
	load_city_people_num();
	initial();
	load_struct_people();//y
	load_I_pro_ini();//y
	generateNetwork_ini();//y
	load_Matrix_w();

	load_write_file();
	return 0;
}
void load_city_people_num()
{
	city_people_num = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	int i, a;
	errno_t err;
	err = fopen_s(&fp, "peoplenum_a.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		fscanf_s(fp, "%d ", &a);
		city_people_num[i] = a;
	}
	fclose(fp);
}
void initial()
{
	adjacentMatrix_w = (float**)malloc(sizeof(float *) * NETWORK_SIZE); //分配指针数组
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_w[i] = (float *)malloc(sizeof(float) * NETWORK_SIZE);//分配每个指针指向的数组
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
	peopleTodes_num = (int**)malloc(sizeof(int *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		peopleTodes_num[i] = (int *)malloc(sizeof(int) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	I_pro_ini = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	P = (double *)malloc(sizeof(double)*NETWORK_SIZE );
	PI = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	rho = (double *)malloc(sizeof(double)*NETWORK_SIZE);
}
void load_struct_people()
{
	int i, j,m;
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
				city[i][j].ini_state = 'S';
			city[i][j].des = 22;
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
void generateNetwork_ini()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
		for (j = i; j < NETWORK_SIZE; j++)
			adjacentMatrix_w[i][j] = adjacentMatrix_w[j][i] = 0;//初始化ER网络的邻接矩阵,全部置为0
	int count = 0;//用以统计网络中无向边的个数
	double probability = 0.0;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = i + 1; j < NETWORK_SIZE; j++)
		{
			probability = rand() / (RAND_MAX + 0.0);//生成一个随机数
			if (probability < PROBABILITY_OF_EAGE)//如果此随机数小于连边概率，则在此（i，j)节点对之间添加一条边，否则不添加边。
			{
				count++;
				adjacentMatrix_w[i][j] = adjacentMatrix_w[j][i] = 1;
			}
		}
	}//重复直到所有的节点对都被选择一次
}
void load_Matrix_w()
{
	int i, j;
	float a;
	errno_t err;
	err = fopen_s(&fp5, "w_a.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			fscanf_s(fp5, "%f ", &a);
			adjacentMatrix_w[i][j] = a;
		}
	}
	fclose(fp5);
}

void load_Matrix_R_fenzi()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] == 0 || rho[j] == 0)
			{
				adjacentMatrix_R[i][j] = 0;
			}
			else
				adjacentMatrix_R[i][j] = pow(adjacentMatrix_w[i][j] * rho[j], la);
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
			if (sum_w[i] == 0)
				adjacentMatrix_R[i][j] = 0;
			else
				adjacentMatrix_R[i][j] = adjacentMatrix_R[i][j] / sum_w[i];//出现误差
		}
	}
}
void load_peopleTodes_num()
{
	int i, j;
	for (i= 0; i< NETWORK_SIZE ; i++)
	{
		for (j= 0; j < NETWORK_SIZE ; j++)
		{
			if (i == j)
			{
				peopleTodes_num[j][i] = (int)((1 - pd)*city_people_num[i] + 0.5);
			}
			if (i != j)
			{
				peopleTodes_num[j][i] = (int)(pd * adjacentMatrix_R[j][i] * city_people_num[j] + 0.5);
			}
		}
	}
}
void load_P()
{
	int i, j;
	double a;
	for (i = 0; i < NETWORK_SIZE ; i++)
	{
		a = 1;
		for (j = 0; j < NETWORK_SIZE ; j++)
		{
			a = (pow(1 - lambda * rho[i], peopleTodes_num[j][i]))*a;
		}
		P[i] = 1 - a;
	}
}
void load_PI()
{
	int i, j;
	double sum;
	for (i = 0; i < NETWORK_SIZE ; i++)
	{
		sum = 0;
		for (j = 0; j < NETWORK_SIZE ; j++)
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
			rho[i] = (1 - mu)*rho[i] + (1 - rho[i])*PI[i];
	}
}
void load_limit_rho()
{
	int t, i,j;
	double pro;
	for (t = 0; t < 100; t++)
	{
		//printf("la:%If pd:%If t:%d \n", la, pd,t);
		load_rho(t);
		/*printf("rho:\n");
		for (i = 0; i < NETWORK_SIZE; i++)
		{
			printf("%If ", rho[i]);
		}
		printf("\n");
		printf("\n");
		printf("W:\n");
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%f ", adjacentMatrix_w[j][i]);
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
				printf("%If ", adjacentMatrix_R[j][i]);
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
			printf("\n");
		}
		printf("\n");*/
		load_peopleTodes_num();
		load_P();
		load_PI();
	}
}

double load_I_pro_zong()
{
	double pro_zong;
	int I_zong = 0, zong = 0;
	int i;
	for (i = 0; i < NETWORK_SIZE ; i++)
	{
		I_zong = (int)(rho[i] * city_people_num[i] + 0.5) + I_zong;
		zong = city_people_num[i] + zong;
	}
	//I_zong = (int)(I_zong + 0.5);
	pro_zong = (double)I_zong / zong;
	return pro_zong;
}
void load_write_file()
{
	double pro;
	int i, j;
	errno_t err;
	err = fopen_s(&fp4, "result1_a.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (pd = 0; pd <= 1;pd=pd+0.01)//11
	{
		for (lambda = 0; lambda <=0.01; lambda=lambda+0.0001)//101
		{
			//for (la = -30; la <= 30; la = la + 0.6)//101
		    //{
			load_limit_rho();
			pro = load_I_pro_zong();
			fprintf_s(fp4, "%lf ", pro);
			//printf("%f %f %f %.2f", pd, la, lambda, pro);
		    //}
			
		   }
		fprintf_s(fp4, "\n");
		//printf("\n");
	}
	fclose(fp4);
}

/*for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%f ", adjacentMatrix[j][i]);
			}
			printf("\n");
		}
		printf("\n");*/

/*
printf("rho:\n");
		for (i = 0; i < NETWORK_SIZE; i++)
		{
			printf("%f ", rho[i]);
		}
		printf("\n");
		printf("\n");
*/

/*
void load_city_people_num()
{
	city_people_num = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	int i, a, num;
	errno_t err;
	err = fopen_s(&fp, "people3.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		num = 0;
		while (1)
		{
			fscanf_s(fp, "%d", &a);
			if (feof(fp))
				break;
			if (a == i)
			{
				num++;
				fseek(fp, 3, 1);
			}
			else
				break;
		}
		city_people_num[i] = num;
	}
	fclose(fp);
}
*/

/*
void load_struct_people()
{
	int i, j;
	errno_t err;
	err = fopen_s(&fp3, "people3.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < city_people_num[i]; j++)
		{
			fscanf_s(fp3, "%d", &city[i][j].sta);
			fscanf_s(fp3, "%c", &city[i][j].ini_state, sizeof(char));
			city[i][j].des = 22;
			city[i][j].fin_state = 'N';
		}
	}
	fclose(fp3);
}
*/