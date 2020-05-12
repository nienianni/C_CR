#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
int NETWORK_SIZE = 11;
int K=6;
double P1 = 0.5;
float ** adjacentMatrix;
float ** adjacentMatrix_help;
//int Nmax = 300;
float delta = 0.4;
//float Alpha = 0.1;
float pd;
float lambda;
float mu = 0.2;
float la = 1;
FILE * fp3;
FILE * fp4;
FILE *fp;
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

void load_city_people_num();
void initial();
void generate_near_Network_ini();
void generate_WS_Network();
void load_Matrix_help();
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
float load_I_pro_zong();
void load_write_file();

int main()
{
	int i, j;
	load_city_people_num();
	initial();
	load_struct_people();//y
	load_I_pro_ini();//y
	generate_near_Network_ini();
	generate_WS_Network();//y
	load_Matrix_help();

	load_write_file();
	return 0;
}
void load_city_people_num()
{
	city_people_num = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	int i, a, num;
	errno_t err;
	err = fopen_s(&fp, "people1.txt", "r");
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
void initial()
{
	adjacentMatrix = (float**)malloc(sizeof(float *) * NETWORK_SIZE); //分配指针数组
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix[i] = (float *)malloc(sizeof(float) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	adjacentMatrix_help = (float**)malloc(sizeof(float *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_help[i] = (float *)malloc(sizeof(float) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	city = (struct people**)malloc(sizeof(struct people *) * NETWORK_SIZE);//有多少个城市,city为指向指针的指针
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		city[i] = (struct people*)malloc(sizeof(struct people)* city_people_num[i]);//叶子城市里有多少个人
	}
	sum_w = (float *)malloc(sizeof(float)*NETWORK_SIZE);
	peopleTodes_num = (int**)malloc(sizeof(int *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		peopleTodes_num[i] = (int *)malloc(sizeof(int) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	I_pro_ini = (float *)malloc(sizeof(float)*NETWORK_SIZE);
	P = (float *)malloc(sizeof(float)*NETWORK_SIZE);
	PI = (float *)malloc(sizeof(float)*NETWORK_SIZE);
	rho = (float *)malloc(sizeof(float)*NETWORK_SIZE);
}
void load_struct_people()
{
	int i, j;
	errno_t err;
	err = fopen_s(&fp3, "people1.txt", "r");
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
		I_pro_ini[m] = (float)I / (I + S);
	}
}
void generate_near_Network_ini()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
		for (j = 0; j < NETWORK_SIZE; j++)
			adjacentMatrix[i][j] = 0;//初始化网络的邻接矩阵,全部置为0
	for (i = 0; i < NETWORK_SIZE; i++)//i代表出发点
	{
		for (j = 1; j <= K / 2; j++)//j代表走几步个邻居
		{
			if (i - j >= 0 && i + j < NETWORK_SIZE)//若朝左朝右走后 都未越界，左右各连一个邻居
			{
				adjacentMatrix[i][i - j] = adjacentMatrix[i][i + j] = 1;
			}
			else if (i - j < 0)//若朝左走后到越下界，左右各连一个邻居
			{
				adjacentMatrix[i][NETWORK_SIZE + i - j] = adjacentMatrix[i][i + j] = 1;
			}
			else if (i + j >= NETWORK_SIZE)//若朝右走后越上界，左右各连一个邻居
			{
				adjacentMatrix[i][i + j - NETWORK_SIZE] = adjacentMatrix[i][i - j] = 1;
			}
		}
	}
}
void generate_WS_Network()
{
	int i, j, *hasEage;
	double isChange = 0.0;//与P比较，决定是否重连
	int re_connectRandomNode;//重连的节点
	hasEage = (int *)malloc(sizeof(int) * NETWORK_SIZE);//用于在选择另外一个节点进行重连的时候，判断是否重边，以及是否随机选择了自身
	int number_changedEage = 0;
	for (i = 0; i < NETWORK_SIZE; i++)//顺时针遍历每个节点
	{
		for (j = 1; j <= K / 2; j++)
		{
			//isChange = (rand() % NETWORK_SIZE) / (double)NETWORK_SIZE;
			isChange = rand() / (RAND_MAX + 0.0);
			if (isChange < P1)//这条边要重连
			{
				while (1)//前提是要有可重连的点。否则会陷入循环
				{
					re_connectRandomNode = (rand() % NETWORK_SIZE);
					if (adjacentMatrix[i][re_connectRandomNode] == 0 && re_connectRandomNode != i)
						break;
				}
				if (i + j < NETWORK_SIZE)
				{
					adjacentMatrix[i][i + j] = adjacentMatrix[i + j][i] = 0;
				}
				else
				{
					adjacentMatrix[i][i + j - NETWORK_SIZE] = adjacentMatrix[i + j - NETWORK_SIZE][i] = 0;
				}
				adjacentMatrix[i][re_connectRandomNode] = adjacentMatrix[re_connectRandomNode][i] = 1;
				number_changedEage++;
			}
			else
			{
				printf("(%d, %d) no change\n", i, i + j);
			}
		}
	}
}
void load_Matrix_help()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
		for (j = 0; j < NETWORK_SIZE; j++)
			adjacentMatrix_help[i][j] = adjacentMatrix[i][j];
}

void load_Matrix_w()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
		for (j = 0; j < NETWORK_SIZE; j++)
			adjacentMatrix[i][j] = adjacentMatrix_help[i][j];
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (i == 0 && adjacentMatrix[i][j] == 1)
			{
				adjacentMatrix[i][j] = 1.0 / (NETWORK_SIZE - 1);
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
void load_Matrix_R_fenzi()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			adjacentMatrix[i][j] = pow(adjacentMatrix[i][j] * rho[j], la);
		}
	}
}
void load_Matrix_R()
{
	int i, j;
	float s;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		s = 0;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			s = adjacentMatrix[i][j] + s;
		}
		sum_w[i] = s;
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (sum_w[i] == 0)
				adjacentMatrix[i][j] = 0;
			else
				adjacentMatrix[i][j] = adjacentMatrix[i][j] / sum_w[i];//出现误差
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
				peopleTodes_num[j][i] = (int)(pd * adjacentMatrix[j][i] * city_people_num[j] + 0.5);
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
	float sum;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		sum = 0;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			sum = (adjacentMatrix[i][j] * P[j]) + sum;
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
	float pro;
	for (t = 0; t < 100; t++)
	{
		load_rho(t);
		load_Matrix_w();
		load_Matrix_R_fenzi();
		load_Matrix_R();
		load_peopleTodes_num();
		load_P();
		load_PI();
	}
}

float load_I_pro_zong()
{
	float pro_zong;
	int I_zong = 0, zong = 0;
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
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
	int i, j;
	errno_t err;
	err = fopen_s(&fp4, "result1_a.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (pd = 0; pd <= 1; pd = pd + 0.01)
	{
		for (lambda = 0; lambda <= 0.003; lambda = lambda + 0.00003)
		{
			load_limit_rho();
			pro = load_I_pro_zong();
			fprintf_s(fp4, "%f ", pro);
			printf("%f  %f  %.6f", pd, lambda, pro);
		}
		fprintf_s(fp4, "\n");
		printf("\n");
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