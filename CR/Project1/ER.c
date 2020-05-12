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
int *wSum_loc;//��¼W����ÿһ�������ֵ��ÿһ�еĺͣ���λ��

float *I_pro_ini;
float *I_pro;
int ** peopleTodes_num;//��¼ÿ�����е�ÿ�����е��������������ڼң�
float *P;//��¼��ɢ������S̬������ÿ����Ⱥ�лᱻ��Ⱦ�ĸ���

FILE * fp1;
FILE * fp2;
int K;//Ҷ����
int Nmax = 300;
float delta = 0.4;//Ҷ��ȥ���ĵĸ���
float Alpha = 0.1;
float pd;//���ŵĸ���
float lambda;
float mu = 0.2;

void initial();//��ʼ������
void generate_ER_Network_ini();//�����ڽӾ���ÿ���ڵ�ɳ���Щ�ط���
void load_struct_people();//���ļ� ����.sta��.ini_state��ʼ��people�ṹ��
void load_Matrix_w();//���ڽӾ����Ϊ��Ȩ��
void load_Matrix_wSum();//����Ȩ�ڽӾ���Ȩֵ��Ϊ�����õ�ĺ�
float load_randnum();//����[0,1]֮��������
void load_I_pro_ini();

void load_I_pro(int t);//���ɳ�ʼʱ ÿ����Ⱥ�ĸ�Ⱦ�˿�ռ��
void load_people_des();
void load_peopleTodes_num();//��һ�������¼ĳ��ȥĳ�ص�����
void load_P();//��һ�������¼ ��ĳ��S̬�ڵ㱻��Ⱦ�ĸ���
void load_people_current_state(int t);
void load_people_fin_state();
void load_I_pro_limit();//100��������ÿ����Ⱥ�ĸ�Ⱦ��

float load_I_pro_zong();//100�������������Ⱦ��
void load_write_file();


int main()
{
	printf("������Ҷ�ӽڵ����");
	scanf_s("%d", &NETWORK_SIZE);
	printf("���������߸���");
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
	adjacentMatrix = (float**)malloc(sizeof(float *) * (NETWORK_SIZE + 1)); //����ָ������
	int i;
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		adjacentMatrix[i] = (float *)malloc(sizeof(float) * (NETWORK_SIZE + 1));//����ÿ��ָ��ָ�������
	}
	city = (struct people**)malloc(sizeof(struct people *) * (NETWORK_SIZE + 1));//�ж��ٸ�����,cityΪָ��ָ���ָ��
	city[0] = (struct people*)malloc(sizeof(struct people)*Nmax);//���ĳ������ж��ٸ��ˣ�ָ��ṹ���ָ��
	for (i = 1; i < NETWORK_SIZE + 1; i++)
	{
		city[i] = (struct people*)malloc(sizeof(struct people)* Alpha*Nmax);//Ҷ�ӳ������ж��ٸ���
	}
	wSum_loc = (int *)malloc(sizeof(int)*(NETWORK_SIZE + 1));
	I_pro_ini = (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));
	I_pro = (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));

	peopleTodes_num = (int**)malloc(sizeof(int *) * (NETWORK_SIZE + 1)); //����ָ������
	for (i = 0; i < NETWORK_SIZE + 1; i++)
	{
		peopleTodes_num[i] = (int *)malloc(sizeof(int) * (NETWORK_SIZE + 1));//����ÿ��ָ��ָ�������
	}
	P = (float *)malloc(sizeof(float)*(NETWORK_SIZE + 1));
}
void generate_ER_Network_ini()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE+1; i++)
		for (j = i; j < NETWORK_SIZE+1; j++)
			adjacentMatrix[i][j] = adjacentMatrix[j][i] = 0;//��ʼ��ER������ڽӾ���,ȫ����Ϊ0
	int count = 0;//����ͳ������������ߵĸ���
	double probability = 0.0;
	for (i = 0; i < NETWORK_SIZE+1; i++)
	{
		for (j = i + 1; j < NETWORK_SIZE+1; j++)
		{
			probability = rand() / (RAND_MAX + 0.0);//����һ�������
			printf("%f ", probability);
			if (probability < PROBABILITY_OF_EAGE)//����������С�����߸��ʣ����ڴˣ�i��j)�ڵ��֮�����һ���ߣ�������ӱߡ�
			{
				count++;
				adjacentMatrix[i][j] = adjacentMatrix[j][i] = 1;
			}
		}
	}//�ظ�ֱ�����еĽڵ�Զ���ѡ��һ��
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
		puts("���ܴ��ļ�");
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
		S = 0;//��¼ÿ����Ⱥ��S����
		I = 0;//��¼ÿ����Ⱥ��I����
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
			if (r_pd < pd)//Ҫ��ȥ
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
			else//���ڼ�
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
		S = 0;//��¼ÿ����Ⱥ��S����
		I = 0;//��¼ÿ����Ⱥ��I����
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

float load_I_pro_zong()//�������I_pro
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
		puts("���ܴ��ļ�");
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