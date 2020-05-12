#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
//10.4 0.1  20.9 0.95  30.4 0.4
int NETWORK_SIZE;
double P1;
int K;
float ** adjacentMatrix;
float ** adjacentMatrix_help;

struct people
{
	int sta;
	char ini_state;
	int des;
	char fin_state;
	char current_state;
}** city;
int *city_people_num;//记录每个城市的初始人口数；
int *wSum_loc;//记录W矩阵每一行最大数值（每一行的和）的位置

float *I_pro_ini;
float *I_pro;
int ** peopleTodes_num;//记录每个城市到每个城市的人数（包括待在家）
float *P;//记录扩散发生后，S态个体在每个种群中会被感染的概率

FILE * fp1;
FILE * fp2;
FILE *fp;
float delta = 0.4;//叶子去中心的概率
//float Alpha = 0.1;
float pd;//出门的概率
float lambda;
float mu = 0.2;
float la = 1;

void load_city_people_num();
void initial();//初始化数组
void load_struct_people();//读文件 读入.sta和.ini_state初始化people结构体
void generate_near_Network_ini();
void generate_WS_Network();
void load_Matrix_help();
void load_I_pro_ini();
void load_I_pro(int t);//生成初始时 每个种群的感染人口占比
void load_Matrix_w();//将邻接矩阵变为带权的
void load_Matrix_R();
void load_Matrix_RSum();//将带权邻接矩阵权值变为：到该点的和
float load_randnum();//生成[0,1]之间的随机数

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
	int i, j;
	printf("请输入节点个数");
	scanf_s("%d", &NETWORK_SIZE);
	do {
		printf("请输入K值");
		scanf_s("%d", &K);
	} while (K % 2 != 0);
	printf("请输入重连概率");
	scanf_s("%lf", &P1);
	load_city_people_num();
	initial();
	load_struct_people();
	load_I_pro_ini();
	generate_near_Network_ini();
	generate_WS_Network();
	load_Matrix_help();

	load_write_file();
	return 0;
}
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
		city[i] = (struct people*)malloc(sizeof(struct people)*city_people_num[i]);//叶子城市里有多少个人
	}
	wSum_loc = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	I_pro_ini = (float *)malloc(sizeof(float)*NETWORK_SIZE);
	I_pro = (float *)malloc(sizeof(float)*NETWORK_SIZE);

	peopleTodes_num = (int**)malloc(sizeof(int *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		peopleTodes_num[i] = (int *)malloc(sizeof(int) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	P = (float *)malloc(sizeof(float)*NETWORK_SIZE);
}
void load_struct_people()
{
	int i, j;
	errno_t err;
	err = fopen_s(&fp1, "people3.txt", "r");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < city_people_num[i]; j++)
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
				if (city[i][j].ini_state == 'S'&&city[i][j].sta == m)
					S++;
				if (city[i][j].ini_state == 'I'&&city[i][j].sta == m)
					I++;
			}
		}
		I_pro_ini[m] = ((float)I / (I + S));
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
float load_randnum()
{
	float Rnum;
	Rnum = rand() / (RAND_MAX + 1.0);
	return Rnum;
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
void load_Matrix_R()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			adjacentMatrix[i][j] = pow(adjacentMatrix[i][j] * I_pro[j], la);
		}
	}
}
void load_Matrix_RSum()
{
	int i, j;
	int loc = 0;
	float pre_sum;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		pre_sum = 0.0;
		for (j = 0; j < NETWORK_SIZE; j++)
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
void load_people_des()
{
	int i, j, m;
	float r_pd;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < city_people_num[i]; j++)
		{
			r_pd = load_randnum();
			if (r_pd < pd)//要出去
			{
				float r_where;
				r_where = (adjacentMatrix[city[i][j].sta][wSum_loc[i]])* (rand() / (RAND_MAX + 0.0));
				for (m = 0; m < NETWORK_SIZE; m++)
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
	int i, j, m, n, Num;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			Num = 0;
			for (m = 0; m < NETWORK_SIZE; m++)
			{
				for (n = 0; n < city_people_num[m]; n++)
				{
					if (city[m][n].sta == i && city[m][n].des == j)
						Num++;
				}
			}
			peopleTodes_num[i][j] = Num;
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
			a = (pow(1 - lambda * I_pro[j], peopleTodes_num[j][i]))*a;
		}
		P[i] = 1 - a;
	}
}
void load_people_current_state(int t)
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < city_people_num[i]; j++)
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
	int i, j;
	float r1, r2;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < city_people_num[i]; j++)
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
void load_I_pro(int t)
{
	int i, j, m, S, I;
	for (m = 0; m < NETWORK_SIZE; m++)
	{
		S = 0;//记录每个种群的S人数
		I = 0;//记录每个种群的I人数
		if (t == 0)
			I_pro[m] = I_pro_ini[m];
		else
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				for (j = 0; j < city_people_num[i]; j++)
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
void load_I_pro_limit()
{
	int i, j, t;
	for (t = 0; t < 100; t++)
	{
		load_I_pro(t);
		load_Matrix_w();
		load_Matrix_R();
		load_Matrix_RSum();
		load_people_des();
		load_peopleTodes_num();
		load_P();
		load_people_current_state(t);
		load_people_fin_state();
	}
}

float load_I_pro_zong()//求总体的I_pro
{
	int I_zong = 0, zong = 0, i;
	float pro_zong;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		I_zong = (int)(I_pro[i] * city_people_num[i] + 0.5) + I_zong;
		zong = city_people_num[i] + zong;
	}
	pro_zong = ((float)I_zong / zong);
	return pro_zong;
}
void load_write_file()
{
	float pro;
	errno_t err;
	err = fopen_s(&fp2, "result3.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (pd = 0; pd <= 1; pd = pd + 0.01)
	{
		for (lambda = 0; lambda <= 0.0027; lambda = lambda + 0.000027)
		{
			load_I_pro_limit();
			pro = load_I_pro_zong();
			fprintf_s(fp2, "%f ", pro);
			printf("%f %f %f", pd, lambda, pro);
		}
		fprintf_s(fp2, "\n");
		printf("\n");
	}
	fclose(fp2);
}