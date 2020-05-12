#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
int NETWORK_SIZE_ini = 30;//用于搭建网络
double PROBABILITY_OF_EAGE = 0.5;
int NETWORK_SIZE = 50;
int e_link_num=15;//每次新加入节点连接的边数
float ** adjacentMatrix_w;
int *sum_degree;//记录每个节点的度数
double ** adjacentMatrix_R;

double pd=0.5;
double lambda;
double mu = 0.2;
double la ;
FILE * fp4;
FILE * fp5;
FILE *fp6;
int *city_people_num;
double *sum_w;
int ** peopleTodes_num;
double*I_pro_ini;
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
void generate_SF_Network_ini();//生成SF网络
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
	generate_SF_Network_ini();
	load_Matrix_w();

	load_write_file();
	return 0;
}
void load_city_people_num()
{
	city_people_num = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	int i, a;
	errno_t err;
	err = fopen_s(&fp4, "peoplenum_a.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		fscanf_s(fp4, "%d ", &a);
		city_people_num[i] = a;
	}
	fclose(fp4);
}
void initial()
{
	adjacentMatrix_w = (float**)malloc(sizeof(float *) * NETWORK_SIZE); //分配指针数组
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_w[i] = (float *)malloc(sizeof(float) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	sum_degree = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	adjacentMatrix_R = (double**)malloc(sizeof(double *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_R[i] = (double *)malloc(sizeof(double) * NETWORK_SIZE);//分配每个指针指向的数组
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
	I_pro_ini = (double*)malloc(sizeof(double)*NETWORK_SIZE);
	P = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	PI = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	rho = (double *)malloc(sizeof(double)*NETWORK_SIZE);
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
	for (i = 0; i < NETWORK_SIZE_ini; i++)
	{
		for (j = i + 1; j < NETWORK_SIZE_ini; j++)
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
void generate_SF_Network_ini()
{
	int i, j, m, n, k, l, sum_e = 0, sum_z = 0;
	double r;
	for (m = 0; m < NETWORK_SIZE; m++)//初始化存放每个节点度数的数组
	{
		sum_degree[m] = 0;
	}
	for (m = NETWORK_SIZE_ini; m < NETWORK_SIZE; m++)//依次加入每个节点
	{
		for (i = 0; i < m; i++)//得到前m个节点的度
		{
			sum_e = 0;
			for (j = 0; j < m; j++)
			{
				sum_e = adjacentMatrix_w[i][j] + sum_e;
			}
			sum_degree[i] = sum_e;
		}
		for (i = 1; i < m; i++)//节点度数和
		{
			sum_degree[i] = sum_degree[i - 1] + sum_degree[i];
		}
		for (l = 0; l < e_link_num; l++)//依次加入节点，并与e_link_num个节点相连
		{
			r = sum_degree[m - 1] * (rand() / (RAND_MAX + 0.0));
			for (k = 0; k < m; k++)
			{
				if (sum_degree[k] >= r)
				{
					adjacentMatrix_w[m][k] = adjacentMatrix_w[k][m] = 1;
					break;
				}
			}
		}
	}
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
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
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
	float a;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		a = 1;
		for (j = 0; j < NETWORK_SIZE; j++)
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
			rho[i] = (1 - mu)*rho[i] + (1 - rho[i])*PI[i];
	}
}
void load_limit_rho()
{
	int t, i, j;
	double pro;
	for (t = 0; t < 100; t++)
	{
		load_rho(t);
		load_Matrix_R_fenzi();
		load_Matrix_R();
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
	for (i = 0; i < NETWORK_SIZE; i++)
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
	err = fopen_s(&fp6, "I_L_la_a.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (la=-5;la<=5;la++)
	{
		for (lambda = 0; lambda <= 0.02; lambda = lambda + 0.0002)
		{
			load_limit_rho();
			pro = load_I_pro_zong();
			fprintf_s(fp6, "%If ", pro);
			//printf("%f  %f  %.6f", pd, lambda, pro);
		}
		fprintf_s(fp6, "\n");
		//printf("\n");
	}
	fclose(fp6);
}