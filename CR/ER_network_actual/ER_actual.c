#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
int NETWORK_SIZE = 10;//叶子数
double PROBABILITY_OF_EAGE = 0.5;
float ** adjacentMatrix;
int Nmax = 300;
float delta = 0.4;
float Alpha = 0.1;
float pd;
float lambda;
float mu = 0.2;
int K;
FILE * fp3;
FILE * fp4;
int *city_people_num;
float *sum_w;
int ** peopleTodes_num;
float *I_pro_ini;
float *P;
float *PI;
float *rho;

struct people
{
	int sta;
	char ini_state;
	int des;
	char fin_state;
}** city;

void initial();
void generateNetwork_ini();
void load_Matrix_w();//将邻接矩阵变为带权的
void load_Matrix_R();
void load_struct_people();
void load_city_people_num();
void load_peopleTodes_num();
void load_I_pro_ini();
void load_P();
void load_PI();
void load_rho(int t);//100个步长后的每个集群的感染占比
void load_limit_rho();
float load_I_pro_zong();
void load_write_file();

int main()
{
	int i, j, N;
	initial();
	generateNetwork_ini();//y
	load_Matrix_w();//y
	load_Matrix_R();
	load_struct_people();//y
	load_city_people_num();//y
	load_I_pro_ini();//y

	load_write_file();
	return 0;
}

void initial()
{
	adjacentMatrix = (float**)malloc(sizeof(float *) * (NETWORK_SIZE + 1)); //分配指针数组
	int i;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		adjacentMatrix[i] = (float *)malloc(sizeof(float) * (NETWORK_SIZE + 1));//分配每个指针指向的数组
	}
	city = (struct people**)malloc(sizeof(struct people *) * (NETWORK_SIZE + 1));//有多少个城市,city为指向指针的指针
	city[0] = (struct people*)malloc(sizeof(struct people)*Nmax);//中心城市里有多少个人，指向结构体的指针
	for (i = 1; i < NETWORK_SIZE + 1; i++)
	{
		city[i] = (struct people*)malloc(sizeof(struct people)* Alpha*Nmax);//叶子城市里有多少个人
	}
	city_people_num = (int *)malloc(sizeof(int)*(NETWORK_SIZE + 1));
	sum_w= (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));
	peopleTodes_num = (int**)malloc(sizeof(int *) * (NETWORK_SIZE + 1)); //分配指针数组
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		peopleTodes_num[i] = (int *)malloc(sizeof(int) * (NETWORK_SIZE + 1));//分配每个指针指向的数组
	}
	I_pro_ini = (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));
	P = (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));
	PI = (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));
	rho = (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));
}
void generateNetwork_ini()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
		for (j = i; j < NETWORK_SIZE + 1; j++)
			adjacentMatrix[i][j] = adjacentMatrix[j][i] = 0;//初始化ER网络的邻接矩阵,全部置为0
	int count = 0;//用以统计网络中无向边的个数
	double probability = 0.0;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		for (j = i + 1; j < NETWORK_SIZE + 1; j++)
		{
			probability = rand() / (RAND_MAX + 0.0);//生成一个随机数
			if (probability < PROBABILITY_OF_EAGE)//如果此随机数小于连边概率，则在此（i，j)节点对之间添加一条边，否则不添加边。
			{
				count++;
				adjacentMatrix[i][j] = adjacentMatrix[j][i] = 1;
			}
		}
	}//重复直到所有的节点对都被选择一次
}
void load_Matrix_w()
{
	int i, j;
	K = NETWORK_SIZE;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		for (j = 0; j < NETWORK_SIZE + 1; j++)
		{
			if (i == 0 && adjacentMatrix[i][j] == 1)
			{
				adjacentMatrix[i][j] = (1.0 / K);
			}
			if (i != 0 && j == 0 && adjacentMatrix[i][j] == 1)
			{
				adjacentMatrix[i][j] = delta;
			}
			if (adjacentMatrix[i][j] == 1)
			{
				adjacentMatrix[i][j] = 1 - delta;
			}
		}
	}
}
void load_Matrix_R()
{
	int i, j;
	float s;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		s = 0;
		for (j = 0; j < NETWORK_SIZE + 1; j++)
		{
			s = adjacentMatrix[i][j] + s;
		}
		sum_w[i] = s;
	}

	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		for (j = 0; j < NETWORK_SIZE + 1; j++)
		{
			 adjacentMatrix[i][j] =adjacentMatrix[i][j]/sum_w[i];//出现误差
		}
	}
}
void load_struct_people()
{
	int i, j, N;
	errno_t err;
	err = fopen_s(&fp3, "people1.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		if (i == 0)
			N = Nmax;
		else
			N = (int)(Alpha * Nmax);
		for (j = 0; j < N; j++)
		{
			fscanf_s(fp3, "%d", &city[i][j].sta);
			fscanf_s(fp3, "%c", &city[i][j].ini_state, sizeof(char));
			city[i][j].des = 22;
			city[i][j].fin_state = 'N';
		}
	}
	fclose(fp3);
}
void load_city_people_num()
{
	int m, i, j, N, num;
	for (m = 0; m < NETWORK_SIZE + 1; m++)
	{
		num = 0;
		for (i = 0; i < NETWORK_SIZE + 1; i++)
		{
			if (i == 0)
				N = Nmax;
			else
				N = (int)(Alpha * Nmax);
			for (j = 0; j < N; j++)
			{
				if (city[i][j].sta == m)
					num++;
			}
		}
		city_people_num[m] = num;
	}
}
void load_I_pro_ini()
{
	int i, j, m, N, S, I;
	for (m = 0; m < NETWORK_SIZE + 1; m++)
	{
		S = 0;//记录每个种群的S人数
		I = 0;//记录每个种群的I人数
		for (i = 0; i < NETWORK_SIZE + 1; i++)
		{
			if (i == 0)
				N = Nmax;
			else
				N = (int)(Alpha * Nmax);
			for (j = 0; j < N; j++)
			{
				if (city[i][j].sta == m && city[i][j].ini_state == 'S')
					S++;
				if (city[i][j].sta == m && city[i][j].ini_state == 'I')
					I++;
			}
		}
		I_pro_ini[m] = (float)I / (I + S);
	}
}

void load_peopleTodes_num()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		for (j = 0; j < NETWORK_SIZE + 1; j++)
		{
			if (i == j)
			{
				peopleTodes_num[j][i] = (int)((1 - pd)*city_people_num[i] + 0.5);
			}
			if (i != j)
			{
				peopleTodes_num[j][i] = (int)(pd * adjacentMatrix[j][i] * city_people_num[j] + 0.5);
			}
		}
	}
}
void load_P()
{
	int i, j;
	float a;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		a = 1;
		for (j = 0; j < NETWORK_SIZE + 1; j++)
		{
			a = (pow(1 - lambda * rho[i], peopleTodes_num[j][i]))*a;
		}
		P[i] = 1 - a;
	}
}
void load_PI()
{
	int i, j;
	float sum;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		sum = 0;
		for (j = 0; j < NETWORK_SIZE + 1; j++)
		{
			sum = (adjacentMatrix[i][j] * P[j]) + sum;
		}
		PI[i] = (1 - pd)*P[i] + pd * sum;
	}
}
void load_rho(int t)
{
	int i;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		if (t == 0)
			rho[i] = I_pro_ini[i];
		else
			rho[i] = (1 - mu)*rho[i] + (1 - rho[i])*PI[i];
	}
}
void load_limit_rho()
{
	int t, i;
	float pro;
	for (t = 0; t < 100; t++)
	{
		load_rho(t);
		load_P();
		load_PI();
	}
}

float load_I_pro_zong()
{
	float pro_zong;
	int I_zong = 0, zong = 0;
	int i;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		I_zong = (int)(rho[i] * city_people_num[i] + 0.5) + I_zong;
		zong = city_people_num[i] + zong;
	}
	//I_zong = (int)(I_zong + 0.5);
	pro_zong = (float)I_zong / zong;
	return pro_zong;
}
void load_write_file()
{
	float pro;
	int i,j;
	errno_t err;
	err = fopen_s(&fp4, "I_L1a.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (pd = 0; pd <= 1.1; pd = pd + 0.1)
	{
		load_peopleTodes_num();
		for (lambda = 0; lambda <= 0.003 ;lambda = lambda + 0.00003)
		{
			load_limit_rho();
			pro = load_I_pro_zong();
			fprintf_s(fp4, "%f ", pro);
			printf("%f %f:%.6f", pd, lambda, pro);
		}
		fprintf_s(fp4, "\n");
		printf("\n");
	}
	fclose(fp4);
}