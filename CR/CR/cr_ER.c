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
int *city_people_num;//��¼ÿ�����еĳ�ʼ�˿�����
int *wSum_loc;//��¼W����ÿһ�������ֵ��ÿһ�еĺͣ���λ��

double *I_pro_ini;
double *I_pro;
int ** peopleTodes_num;//��¼ÿ�����е�ÿ�����е��������������ڼң�
double *P;//��¼��ɢ������S̬������ÿ����Ⱥ�лᱻ��Ⱦ�ĸ���

//FILE * fp1;
FILE * fp2;
//FILE *fp;
FILE *fp3;
FILE *fp4;
double pd;//���ŵĸ���
double lambda ;
double mu=0.2;
double la=-2;

int rand_range_i(int min, int max);
double normal(int x, int miu, double sigma);
int rand_all_peoplenum(int miu, double sigma, int min, int max);
void load_city_people_num();
void initial();//��ʼ������
void load_struct_people();//���ļ� ����.sta��.ini_state��ʼ��people�ṹ��
void generate_ER_Network_ini();//�����ڽӾ���ÿ���ڵ�ɳ���Щ�ط���
int rand_range_f(int min, int max);
double power_low(int x, double c, double a);
int rand_all_w(double c, double a, int min, int max);
void load_Matrix_w();//���ڽӾ����Ϊ��Ȩ��
void load_I_pro_ini();
void load_I_pro(int t);
void load_Matrix_R();
void load_Matrix_RSum();//����Ȩ�ڽӾ���Ȩֵ��Ϊ�����õ�ĺ�
double load_randnum();//����[0,1]֮��������

void load_people_des();
void load_peopleTodes_num();//��һ�������¼ĳ��ȥĳ�ص�����
void load_P();//��һ�������¼ ��ĳ��S̬�ڵ㱻��Ⱦ�ĸ���
void load_people_current_state(int t);
void load_people_fin_state();
void load_I_pro_limit();//100��������ÿ����Ⱥ�ĸ�Ⱦ��

double load_I_pro_zong();//100�������������Ⱦ��
void load_write_file();


int main()
{
	int i,j;
	printf("������ڵ����");
	scanf_s("%d", &NETWORK_SIZE);
	printf("���������߸���");
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
int rand_range_i(int min, int max)//����max��min��һ���������
{
	return rand() % (max - min + 1) + min;
}
double normal(int x, int miu, double sigma)//���ɷ���f(x)���ϵĺ���
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
		puts("���ܴ��ļ�");
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
	adjacentMatrix_w = (float**)malloc(sizeof(float *) * NETWORK_SIZE); //����ָ������
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_w[i] = (float *)malloc(sizeof(float) * NETWORK_SIZE);//����ÿ��ָ��ָ�������
	}
	/*adjacentMatrix_help = (float**)malloc(sizeof(float *) * NETWORK_SIZE); //����ָ������
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_help[i] = (float *)malloc(sizeof(float) * NETWORK_SIZE);//����ÿ��ָ��ָ�������
	}*/
	adjacentMatrix_R= (double**)malloc(sizeof(double *) * NETWORK_SIZE); //����ָ������
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_R[i] = (double *)malloc(sizeof(double) * NETWORK_SIZE);//����ÿ��ָ��ָ�������
	}
	city = (struct people**)malloc(sizeof(struct people *) * NETWORK_SIZE);//�ж��ٸ�����,cityΪָ��ָ���ָ��
	for (i = 0; i < NETWORK_SIZE ; i++)
	{
		city[i] = (struct people*)malloc(sizeof(struct people)*city_people_num[i]);//Ҷ�ӳ������ж��ٸ���
	}
	wSum_loc = (int *)malloc(sizeof(int)*NETWORK_SIZE );
	I_pro_ini = (double *)malloc(sizeof(double)*NETWORK_SIZE );
	I_pro = (double *)malloc(sizeof(double)*NETWORK_SIZE );

	peopleTodes_num = (int**)malloc(sizeof(int *) * NETWORK_SIZE ); //����ָ������
	for (i = 0; i < NETWORK_SIZE ; i++)
	{
		peopleTodes_num[i] = (int *)malloc(sizeof(int) * NETWORK_SIZE );//����ÿ��ָ��ָ�������
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
		S = 0;//��¼ÿ����Ⱥ��S����
		I = 0;//��¼ÿ����Ⱥ��I����
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
			adjacentMatrix_w[i][j] = adjacentMatrix_w[j][i] = 0;//��ʼ��ER������ڽӾ���,ȫ����Ϊ0
	int count = 0;//����ͳ������������ߵĸ���
	double probability = 0.0;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = i + 1; j < NETWORK_SIZE; j++)
		{
			probability = rand() / (RAND_MAX + 0.0);//����һ�������
			if (probability < PROBABILITY_OF_EAGE)//����������С�����߸��ʣ����ڴˣ�i��j)�ڵ��֮�����һ���ߣ�������ӱߡ�
			{
				count++;
				adjacentMatrix_w[i][j] = adjacentMatrix_w[j][i] = 1;
			}
		}
	}//�ظ�ֱ�����еĽڵ�Զ���ѡ��һ��
}
int rand_range_f(int min, int max)//����max��min��һ�����С��
{
	return rand() % (max - min + 1) + min;
	//return min + (max - min)*rand() / (RAND_MAX + 0.0);
}
double power_low(int x, double c, double a)//���ɷ���f(x)���ϵĺ���
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
		puts("���ܴ��ļ�");
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
			if (r_pd < pd)//Ҫ��ȥ
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
			else//���ڼ�
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
		S = 0;//��¼ÿ����Ⱥ��S����
		I = 0;//��¼ÿ����Ⱥ��I����
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

double load_I_pro_zong()//�������I_pro
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
		puts("���ܴ��ļ�");
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
		puts("���ܴ��ļ�");
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
		puts("���ܴ��ļ�");
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