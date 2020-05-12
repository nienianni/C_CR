#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define pi 3.1415926

int NETWORK_SIZE;
double PROBABILITY_OF_EAGE;
float ** adjacentMatrix_w;
double** adjacentMatrix_R;

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

double *I_pro_ini;
double *I_pro;
int ** peopleTodes_num;//记录每个城市到每个城市的人数（包括待在家）
double *P;//记录扩散发生后，S态个体在每个种群中会被感染的概率

//FILE * fp1;
FILE * fp2;
//FILE *fp;
FILE *fp3;
FILE *fp4;
double pd;//出门的概率
double lambda ;
double mu=0.2;
double la=-2;

int rand_range_i(int min, int max);
double normal(int x, int miu, double sigma);
int rand_all_peoplenum(int miu, double sigma, int min, int max);
void load_city_people_num();
void initial();//初始化数组
void load_struct_people();//读文件 读入.sta和.ini_state初始化people结构体
void generate_ER_Network_ini();//生成邻接矩阵，每个节点可朝哪些地方走
int rand_range_f(int min, int max);
double power_low(int x, double c, double a);
int rand_all_w(double c, double a, int min, int max);
void load_Matrix_w();//将邻接矩阵变为带权的
void load_I_pro_ini();
void load_I_pro(int t);
void load_Matrix_R();
void load_Matrix_RSum();//将带权邻接矩阵权值变为：到该点的和
double load_randnum();//生成[0,1]之间的随机数

void load_people_des();
void load_peopleTodes_num();//用一个矩阵记录某地去某地的人数
void load_P();//用一个数组记录 在某地S态节点被感染的概率
void load_people_current_state(int t);
void load_people_fin_state();
void load_I_pro_limit();//100个步长后，每个集群的感染率

double load_I_pro_zong();//100个步长后总体感染率
void load_write_file();


int main()
{
	int i,j;
	printf("请输入节点个数");
	scanf_s("%d", &NETWORK_SIZE);
	printf("请输入连边概率");
	scanf_s("%lf", &PROBABILITY_OF_EAGE);
	//srand((unsigned)time(NULL));
	load_city_people_num();
	initial();
	load_struct_people();
	load_I_pro_ini();
	generate_ER_Network_ini();
	load_Matrix_w();

	load_write_file();
	return 0;
}
int rand_range_i(int min, int max)//生成max到min的一个随机整数
{
	return rand() % (max - min + 1) + min;
}
double normal(int x, int miu, double sigma)//生成符合f(x)线上的函数
{
	return 1.0 / sqrt(2 * pi) / sigma * exp(-1 * (x - miu)*(x - miu) / (2 * sigma*sigma));
}
int rand_all_peoplenum(int miu, double sigma, int min, int max)
{
	int x, dScope;
	double y;
	do {
		x = rand_range_i(min, max);
		y = normal(x, miu, sigma);
		dScope = rand_range_i(0, (int)normal(miu, miu, sigma));
	} while (dScope - y > 0);
	return x;
}
void load_city_people_num()
{
	int i;
	errno_t err;
	err = fopen_s(&fp3, "peoplenum.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	city_people_num = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		city_people_num[i] = rand_all_peoplenum(100, 0.2, 1, 200);
		fprintf_s(fp3, "%d ", city_people_num[i]);
	}
	fclose(fp3);
}
void initial()
{
	adjacentMatrix_w = (float**)malloc(sizeof(float *) * NETWORK_SIZE); //分配指针数组
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_w[i] = (float *)malloc(sizeof(float) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	/*adjacentMatrix_help = (float**)malloc(sizeof(float *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_help[i] = (float *)malloc(sizeof(float) * NETWORK_SIZE);//分配每个指针指向的数组
	}*/
	adjacentMatrix_R= (double**)malloc(sizeof(double *) * NETWORK_SIZE); //分配指针数组
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_R[i] = (double *)malloc(sizeof(double) * NETWORK_SIZE);//分配每个指针指向的数组
	}
	city = (struct people**)malloc(sizeof(struct people *) * NETWORK_SIZE);//有多少个城市,city为指向指针的指针
	for (i = 0; i < NETWORK_SIZE ; i++)
	{
		city[i] = (struct people*)malloc(sizeof(struct people)*city_people_num[i]);//叶子城市里有多少个人
	}
	wSum_loc = (int *)malloc(sizeof(int)*NETWORK_SIZE );
	I_pro_ini = (double *)malloc(sizeof(double)*NETWORK_SIZE );
	I_pro = (double *)malloc(sizeof(double)*NETWORK_SIZE );

	peopleTodes_num = (int**)malloc(sizeof(int *) * NETWORK_SIZE ); //分配指针数组
	for (i = 0; i < NETWORK_SIZE ; i++)
	{
		peopleTodes_num[i] = (int *)malloc(sizeof(int) * NETWORK_SIZE );//分配每个指针指向的数组
	}
	P = (double *)malloc(sizeof(double)*NETWORK_SIZE );
}
void load_struct_people()
{
	int i, j,m;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j <city_people_num[i]; j++)
		{
			city[i][j].sta = i;
			if (i == 0&&j<10)
			{
				city[i][j].ini_state = 'I';
			}
			else
				city[i][j].ini_state = 'S';
			city[i][j].des = 22;
			city[i][j].fin_state = 'N';
			city[i][j].current_state = 'N';
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
				if (city[i][j].ini_state == 'S'&&city[i][j].sta == m)
					S++;
				if (city[i][j].ini_state == 'I'&&city[i][j].sta == m)
					I++;
			}
		}
		I_pro_ini[m] = ((double)I / (I + S));
	}
}
void generate_ER_Network_ini()
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
int rand_range_f(int min, int max)//生成max到min的一个随机小数
{
	return rand() % (max - min + 1) + min;
	//return min + (max - min)*rand() / (RAND_MAX + 0.0);
}
double power_low(int x, double c, double a)//生成符合f(x)线上的函数
{
	return c / pow(x, a);
}
int rand_all_w(double c, double a, int min, int max)
{
	int x, dScope;
	double y;
	do {
		x = rand_range_f(min, max);
		y = power_low(x, c, a);
		dScope = rand_range_f(0, power_low(min, c, a));
	} while (dScope - y > 0);
	return x;
}
void load_Matrix_w()
{
	int i, j;
	errno_t err;
	err = fopen_s(&fp4, "w.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] == 1)
			{
				adjacentMatrix_w[i][j] = rand_all_w(1, 2, 1, 10);
			}
		}
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			fprintf_s(fp4, "%f ", adjacentMatrix_w[i][j]);
		}
	}
	fclose(fp4);
}
double load_randnum()
{
	double Rnum;
	Rnum = rand() / (RAND_MAX + 0.0);
	return Rnum;
}

void load_Matrix_R()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] == 0 || I_pro[j] == 0)
			{
				adjacentMatrix_R[i][j] = 0;
			}
			else
			adjacentMatrix_R[i][j] = pow(adjacentMatrix_w[i][j] * I_pro[j], la);
		}
	}
}
void load_Matrix_RSum()
{
	int i, j;
	int loc = 0;
	double pre_sum;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		pre_sum = 0.0;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			pre_sum = adjacentMatrix_R[i][j] + pre_sum;
			if (adjacentMatrix_R[i][j] != 0)
			{
				adjacentMatrix_R[i][j] = pre_sum;
				loc = j;
			}
		}
		wSum_loc[i] = loc;
	}
}
void load_people_des()
{
	int i, j,m;
	double r_pd;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < city_people_num[i]; j++)
		{
			r_pd = load_randnum();
			if (r_pd < pd)//要出去
			{
				double r_where;
				r_where = (adjacentMatrix_R[city[i][j].sta][wSum_loc[i]])* (rand() / (RAND_MAX + 0.0));
				for (m = 0; m < NETWORK_SIZE ; m++)
				{
					if (adjacentMatrix_R[city[i][j].sta][m] >= r_where)
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
	int i, j,m, n, Num;
	for (i = 0; i < NETWORK_SIZE ; i++)
	{
		for (j = 0; j < NETWORK_SIZE ; j++)
		{
			Num = 0;
			for (m = 0; m < NETWORK_SIZE ; m++)
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
	double a;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		a = 1;
		for (j = 0; j < NETWORK_SIZE ; j++)
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
	double r1, r2;
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
			I_pro[m] = ((double)I / (I + S));
		}
	}
}
void load_I_pro_limit()
{
	int i,j,t;
	for (t = 0; t < 100; t++)
	{
		//printf("la:%f pd:%f t:%d \n", la, pd, t);
		load_I_pro(t);
		/*printf("rho:\n");
		for (i = 0; i < NETWORK_SIZE; i++)
		{
			printf("%f ", I_pro[i]);
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
		load_Matrix_R();
		/*printf("R_fenzi:\n");
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%f ", adjacentMatrix_R[j][i]);
			}
			printf("\n");
		}
		printf("\n");*/
		load_Matrix_RSum();
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
		load_people_des();
		load_peopleTodes_num();
		load_P();
		load_people_current_state(t);
		load_people_fin_state();
	}
}

double load_I_pro_zong()//求总体的I_pro
{
	int I_zong = 0, zong = 0, i;
	double pro_zong;
	for (i = 0; i < NETWORK_SIZE ; i++)
	{
		I_zong = (int)(I_pro[i] * city_people_num[i] + 0.5) + I_zong;
		zong = city_people_num[i] + zong;
	}
	pro_zong = ((double)I_zong / zong);
	return pro_zong;
}
void load_write_file()
{
	double pro;
	errno_t err;
	err = fopen_s(&fp2, "result2.txt", "w");
	if (err != 0)
	{
		puts("不能打开文件");
	}
	for (pd = 0; pd <= 1;pd=pd+0.01)//11
	{
		for (lambda = 0; lambda <= 0.01; lambda = lambda + 0.0001)//101
		{
		//for (la = -30; la <= 30; la = la + 0.6)//101
		//{
			load_I_pro_limit();
			pro = load_I_pro_zong();
			fprintf_s(fp2, "%lf ", pro);
			//printf("%f %f %f %.2f", pd,la,lambda,pro);
		//}
	    }
		fprintf_s(fp2, "\n");
		//printf("\n");
    }
	fclose(fp2);
}

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
		for (j = 0; j <city_people_num[i]; j++)
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
*/