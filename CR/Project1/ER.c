#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
//10.4 0.1  20.9 0.95  30.4 0.4
int NETWORK_SIZE;
double PROBABILITY_OF_EAGE;
float ** adjacentMatrix;
struct people
{
	int sta;
	char ini_state;
	int des;
	char fin_state;
	char current_state;
}** city;
int *wSum_loc;//记录W矩阵每一行最大数值（每一行的和）的位置

float *I_pro_ini;
float *I_pro;
int ** peopleTodes_num;//记录每个城市到每个城市的人数（包括待在家）
float *P;//记录扩散发生后，S态个体在每个种群中会被感染的概率

FILE * fp1;
FILE * fp2;
int K;//叶子数
int Nmax = 300;
float delta = 0.4;//叶子去中心的概率
float Alpha = 0.1;
float pd;//出门的概率
float lambda;
float mu = 0.2;

void initial();//初始化数组
void generate_ER_Network_ini();//生成邻接矩阵，每个节点可朝哪些地方走
void load_struct_people();//读文件 读入.sta和.ini_state初始化people结构体
void load_Matrix_w();//将邻接矩阵变为带权的
void load_Matrix_wSum();//将带权邻接矩阵权值变为：到该点的和
float load_randnum();//生成[0,1]之间的随机数
void load_I_pro_ini();

void load_I_pro(int t);//生成初始时 每个种群的感染人口占比
void load_people_des();
void load_peopleTodes_num();//用一个矩阵记录某地去某地的人数
void load_P();//用一个数组记录 在某地S态节点被感染的概率
void load_people_current_state(int t);
void load_people_fin_state();
void load_I_pro_limit();//100个步长后，每个集群的感染率

float load_I_pro_zong();//100个步长后总体感染率
void load_write_file();


int main()
{
	printf("请输入叶子节点个数");
	scanf_s("%d", &NETWORK_SIZE);
	printf("请输入连边概率");
	scanf_s("%lf", &PROBABILITY_OF_EAGE);
	initial();
	generate_ER_Network_ini();
	load_Matrix_w();
	load_Matrix_wSum();
	load_struct_people();
	load_I_pro_ini();

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
	wSum_loc = (int *)malloc(sizeof(int)*(NETWORK_SIZE + 1));
	I_pro_ini = (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));
	I_pro = (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));

	peopleTodes_num = (int**)malloc(sizeof(int *) * (NETWORK_SIZE + 1)); //分配指针数组
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		peopleTodes_num[i] = (int *)malloc(sizeof(int) * (NETWORK_SIZE + 1));//分配每个指针指向的数组
	}
	P = (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));
}
void generate_ER_Network_ini()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE+1; i++)
		for (j = i; j < NETWORK_SIZE+1; j++)
			adjacentMatrix[i][j] = adjacentMatrix[j][i] = 0;//初始化ER网络的邻接矩阵,全部置为0
	int count = 0;//用以统计网络中无向边的个数
	double probability = 0.0;
	for (i = 0; i < NETWORK_SIZE+1; i++)
	{
		for (j = i + 1; j < NETWORK_SIZE+1; j++)
		{
			probability = rand() / (RAND_MAX + 0.0);//生成一个随机数
			printf("%f ", probability);
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
void load_Matrix_wSum()
{
	int i, j, loc;
	float pre_sum;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		pre_sum = 0.0;
		for (j = 0; j < NETWORK_SIZE + 1; j++)
		{
			pre_sum = adjacentMatrix[i][j] + pre_sum;
			if (adjacentMatrix[i][j] != 0)
			{
				adjacentMatrix[i][j] = pre_sum;
				loc = j;
			}
		}
		wSum_loc[i] = loc;
	}
}
void load_struct_people()
{
	int i, j, N;
	errno_t err;
	err = fopen_s(&fp1, "people1.txt", "r");
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
			fscanf_s(fp1, "%d", &city[i][j].sta);
			fscanf_s(fp1, "%c", &city[i][j].ini_state, sizeof(char));
			city[i][j].des = 22;
			city[i][j].fin_state = 'N';
			city[i][j].current_state = 'N';
		}
	}
	fclose(fp1);
}
float load_randnum()
{
	float Rnum;
	Rnum = rand() / (RAND_MAX + 1.0);
	return Rnum;
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
				if (city[i][j].ini_state == 'S'&&city[i][j].sta == m)
					S++;
				if (city[i][j].ini_state == 'I'&&city[i][j].sta == m)
					I++;
			}
		}
		I_pro_ini[m] = ((float)I / (I + S));
	}
}

void load_people_des()
{
	int i, j, N, m;
	float r_pd;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		if (i == 0)
			N = Nmax;
		else
			N = (int)(Alpha * Nmax);
		for (j = 0; j < N; j++)
		{
			r_pd = load_randnum();
			if (r_pd < pd)//要出去
			{
				float r_where;
				r_where = (adjacentMatrix[city[i][j].sta][wSum_loc[i]])* (rand() / (RAND_MAX + 0.0));
				for (m = 0; m < NETWORK_SIZE + 1; m++)
				{
					if (adjacentMatrix[city[i][j].sta][m] >= r_where)
					{
						city[i][j].des = m;
						break;
					}
				}
			}
			else//待在家
			{
				city[i][j].des = city[i][j].sta;
			}
		}
	}
}
void load_peopleTodes_num()
{
	int i, j, N, m, n, Num;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		for (j = 0; j < NETWORK_SIZE + 1; j++)
		{
			Num = 0;
			for (m = 0; m < NETWORK_SIZE + 1; m++)
			{
				if (m == 0)
					N = Nmax;
				else
					N = (int)(Alpha * Nmax);
				for (n = 0; n < N; n++)
				{
					if (city[m][n].sta == i && city[m][n].des == j)
						Num++;
				}
			}
			peopleTodes_num[i][j] = Num;
		}
	}
}
void load_I_pro(int t)
{
	int i, j, m, N, S, I;
	for (m = 0; m < NETWORK_SIZE + 1; m++)
	{
		S = 0;//记录每个种群的S人数
		I = 0;//记录每个种群的I人数
		if (t == 0)
			I_pro[m] = I_pro_ini[m];
		else
		{
			for (i = 0; i < NETWORK_SIZE + 1; i++)
			{
				if (i == 0)
					N = Nmax;
				else
					N = (int)(Alpha * Nmax);
				for (j = 0; j < N; j++)
				{
					if (city[i][j].fin_state == 'S'&&city[i][j].sta == m)
						S++;
					if (city[i][j].fin_state == 'I'&&city[i][j].sta == m)
						I++;
				}
			}
			I_pro[m] = ((float)I / (I + S));
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
		for (j = 0; j < NETWORK_SIZE+1; j++)
		{
			a = (pow(1 - lambda * I_pro[j], peopleTodes_num[j][i]))*a;
		}
		P[i] = 1 - a;
	}
}
void load_people_current_state(int t)
{
	int i, j, N;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		if (i == 0)
			N = Nmax;
		else
			N = (int)(Alpha * Nmax);
		for (j = 0; j < N; j++)
		{
			if (t == 0)
				city[i][j].current_state = city[i][j].ini_state;
			else
				city[i][j].current_state = city[i][j].fin_state;
		}
	}
}
void load_people_fin_state()
{
	int i, j, N;
	float r1, r2;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		if (i == 0)
			N = Nmax;
		else
			N = (int)(Alpha * Nmax);
		for (j = 0; j < N; j++)
		{
			if (city[i][j].current_state == 'S')
			{
				r1 = load_randnum();
				if (r1 < P[city[i][j].des])
				{
					city[i][j].fin_state = 'I';
				}
				else
					city[i][j].fin_state = city[i][j].current_state;
			}
			if (city[i][j].current_state == 'I')
			{
				r2 = load_randnum();
				if (r2 < mu)
					city[i][j].fin_state = 'S';
				else
					city[i][j].fin_state = city[i][j].current_state;
			}
		}
	}
}
void load_I_pro_limit()
{
	int t, i;
	for (t = 0; t < 100; t++)
	{
		load_people_des();
		load_peopleTodes_num();
		load_I_pro(t);
		load_P();
		load_people_current_state(t);
		load_people_fin_state();
	}
}

float load_I_pro_zong()//求总体的I_pro
{
	int I_zong = 0, zong = 0, i, N;
	float pro_zong;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		if (i == 0)
			N = Nmax;
		else
			N = (int)(Alpha * Nmax);
		I_zong = (int)(I_pro[i] * N + 0.5) + I_zong;
		zong = N + zong;
	}
	pro_zong = ((float)I_zong / zong);
	return pro_zong;
}
void load_write_file()
{
	float pro;
	errno_t err;
	err = fopen_s(&fp2, "I_L1.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (pd = 0; pd <= 1.1; pd = pd + 0.1)
	{
		for (lambda = 0; lambda <= 0.003; lambda = lambda + 0.00003)
		{
			load_I_pro_limit();
			pro = load_I_pro_zong();
			fprintf_s(fp2, "%f ", pro);
			printf("%f %f:%f", pd, lambda, pro);
		}
		fprintf_s(fp2, "\n");
		printf("\n");
	}
	fclose(fp2);
	}