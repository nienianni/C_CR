#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define pi 3.1415926

int NETWORK_SIZE;
double PROBABILITY_OF_EAGE;
int ** adjacentMatrix_w;
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

int *I_people_num;
double *I_pro_ini;
double *I_pro;
double *I_pro_t;
double *sum_w;
int ** peopleTodes_num;//��¼ÿ�����е�ÿ�����е��������������ڼң�
double *P;//��¼��ɢ������S̬������ÿ����Ⱥ�лᱻ��Ⱦ�ĸ���


FILE * fp1;
FILE * fp2;
FILE *fp5;
FILE *fp3;
FILE *fp4;
double pd;//���ŵĸ���
double lambda ;
double mu=0.2;
double la;
int tmax = 2000;//�ܴ���
int t_test = 1000;//ʵ�����

void load_city_people_num();
void initial();//��ʼ������
void load_struct_people();//���ļ� ����.sta��.ini_state��ʼ��people�ṹ��
void generate_ER_Network_ini();//�����ڽӾ���ÿ���ڵ�ɳ���Щ�ط���
void load_Matrix_w();//���ڽӾ����Ϊ��Ȩ��
void load_I_pro_ini();
void load_I_pro(int t);
void load_Matrix_R_fenzi();
void load_Matrix_R();
void load_Matrix_RSum();//����Ȩ�ڽӾ���Ȩֵ��Ϊ�����õ�ĺ�
double load_randnum();//����[0,1]֮��������

void load_people_des();
void load_people_current_state(int t);
void load_I_people_num();
void load_people_fin_state();
void load_I_pro_t(int t);
void load_I_pro_limit();//100��������ÿ����Ⱥ�ĸ�Ⱦ��

double load_I_pro_zong();//100�������������Ⱦ��
double load_I_sus();
void load_write_file();

int main()
{
	int i,j;
	printf("������ڵ����");
	scanf_s("%d", &NETWORK_SIZE);
	printf("���������߸���");
	scanf_s("%lf", &PROBABILITY_OF_EAGE);
	srand((unsigned)time(NULL));
	load_city_people_num();
	initial();
	load_struct_people();
	load_I_pro_ini();
	//generate_ER_Network_ini();
	load_Matrix_w();
	
	load_write_file();
	return 0;
}
void load_city_people_num()
{
	city_people_num = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	int i, a;
	errno_t err;
	err = fopen_s(&fp1, "peoplenum.txt", "r");
	if (err != 0)
	{
		puts("���ܴ��ļ�");
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
	adjacentMatrix_w = (int**)malloc(sizeof(int *) * NETWORK_SIZE); //����ָ������
	int i;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		adjacentMatrix_w[i] = (int *)malloc(sizeof(int) * NETWORK_SIZE);//����ÿ��ָ��ָ�������
	}
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
	sum_w = (double *)malloc(sizeof(double)*NETWORK_SIZE);
	wSum_loc = (int *)malloc(sizeof(int)*NETWORK_SIZE );
	I_people_num = (int *)malloc(sizeof(int)*NETWORK_SIZE);
	I_pro_ini = (double *)malloc(sizeof(double)*NETWORK_SIZE );
	I_pro = (double *)malloc(sizeof(double)*NETWORK_SIZE );
	I_pro_t = (double *)malloc(sizeof(double)*t_test);

	peopleTodes_num = (int**)malloc(sizeof(int *) * NETWORK_SIZE ); //����ָ������
	for (i = 0; i < NETWORK_SIZE ; i++)
	{
		peopleTodes_num[i] = (int *)malloc(sizeof(int) * NETWORK_SIZE );//����ÿ��ָ��ָ�������
	}
	P = (double *)malloc(sizeof(double)*NETWORK_SIZE );
}
void load_struct_people()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j <city_people_num[i]; j++)
		{
			city[i][j].sta = i;
			if (j< (int)(city_people_num[i] * 0.01 + 0.5))
			{
				city[i][j].ini_state = 'I';
			}
			else
			{
				city[i][j].ini_state = 'S';
			}
			city[i][j].des = 100;
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
/*void generate_ER_Network_ini()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
		for (j = i; j < NETWORK_SIZE; j++)
			adjacentMatrix_w[i][j] = adjacentMatrix_w[j][i] = 0;//��ʼ��ER������ڽӾ���,ȫ����Ϊ0
	int count = 0;//����ͳ������������ߵĸ���
	double probability = 0;
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
	printf("%d\n", count*2);
}
void load_Matrix_w()
{
	int i, j;
	errno_t err;
	errno_t err1;
	int b;
	err = fopen_s(&fp2, "w.txt", "r");
	err1 = fopen_s(&fp3, "w_a.txt", "w");
	if (err != 0)
	{
		puts("���ܴ��ļ�");
	}
	if (err1 != 0)
	{
		puts("���ܴ��ļ�w_a");
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] == 1)
			{
				fscanf_s(fp2, "%d ", &b);
				adjacentMatrix_w[i][j] = b;
			}
		}
	}
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			fprintf_s(fp3, "%d ", adjacentMatrix_w[i][j]);
		}
	}
	fclose(fp2);
	fclose(fp3);
}*/
void load_Matrix_w()
{
	int i, j;
	int a;
	errno_t err;
	err = fopen_s(&fp2, "w_a.txt", "r");
	if (err != 0)
	{
		puts("���ܴ��ļ�");
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
}
double load_randnum()
{
	double Rnum;
	Rnum = rand() / (RAND_MAX + 0.0);
	return Rnum;
}

void load_Matrix_R_fenzi()
{
	int i, j;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			if (adjacentMatrix_w[i][j] == 0 || I_pro[j] ==0)
			{
				adjacentMatrix_R[i][j] = 0;
			}
			else
			{
				adjacentMatrix_R[i][j] =pow(adjacentMatrix_w[i][j] * I_pro[j], la);
			}
		}
	}
}
void load_Matrix_RSum()
{
	int i, j;
	int loc=0 ;
	double pre_sum;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		pre_sum = 0;
		for (j = 0; j < NETWORK_SIZE; j++)
		{
			pre_sum = adjacentMatrix_R[i][j] + pre_sum;
			if (adjacentMatrix_R[i][j] != 0)
			{
				adjacentMatrix_R[i][j] = pre_sum;
				loc = j;
			}
		}
		wSum_loc[i] = loc;//��¼��һ�������ֵ����һ��
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
				r_where = (adjacentMatrix_R[city[i][j].sta][wSum_loc[city[i][j].sta]])* (rand() / (RAND_MAX + 0.0));
				//printf("r_where=%f\n", r_where);
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
void load_I_people_num()//��¼��ɢ��ÿ����Ⱥ�еĸ�Ⱦ����
{
	int i, j,m,num;
	for (m = 0; m < NETWORK_SIZE; m++)
	{
		num = 0;
		for (i = 0; i < NETWORK_SIZE; i++)
		{
			for (j = 0; j < city_people_num[i]; j++)
			{
				if (city[i][j].des == m && city[i][j].current_state == 'I')
				{
					num++;
				}
			}
		}
		I_people_num[m] = num;
	}
}
void load_people_fin_state()
{
	int i, j,m,n,A;
	double r1, r2;
	int I;//��¼�ü�Ⱥ����I̬������
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		for (j = 0; j < city_people_num[i]; j++)
		{
			if (city[i][j].current_state == 'S')
			{
				A = 0;
				for (m = 0; m < I_people_num[city[i][j].des]; m++)
				{
					r1 = load_randnum();
					if (r1 <lambda)
					{
						A = 2;
						city[i][j].fin_state = 'I';
						//break;
					}
				}
				if (A == 0)
				{
					city[i][j].fin_state = 'S';
				}
			}
			if (city[i][j].current_state == 'I')
			{
				r2 = load_randnum();
				if (r2 < mu)
					city[i][j].fin_state = 'S';
				else
					city[i][j].fin_state = 'I';
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
			if (I_pro[m] == 0)
			{
				printf("%d:lalala\n ",m);
			}
		}
	}
}
void load_I_pro_t(int t)
{
	int i;
	double sum = 0;
	for (i = 0; i < NETWORK_SIZE; i++)
	{
		sum = sum + I_pro[i];
	}
	I_pro_t[t-(tmax-t_test)] =sum / (double)NETWORK_SIZE;
}
void load_I_pro_limit()
{
	int i,j,t;
	for (t = 0; t < tmax; t++)
	{
		//printf("la:%f pd:%f t:%d \n", la, pd, t);
		load_I_pro(t);
		printf("rho:\n");
		for (i = 0; i < NETWORK_SIZE; i++)
		{
			printf("%f ", I_pro[i]);
		}
		printf("\n");
		printf("\n");
		if (t >= (tmax-t_test))
		{
			load_I_pro_t(t);
		}
		/*printf("I_pro_t:\n");
		for (i = 0; i < t_test; i++)
		{
			printf("%f ", I_pro_t[i]);
		}
		printf("\n");
		printf("\n");*/
		/*printf("W:\n");
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
		printf("\n");
		*/
		load_people_des();
		/*for (j = 0; j < NETWORK_SIZE; j++)
		{
			for (i = 0; i < NETWORK_SIZE; i++)
			{
				printf("%d ", city[j][i].des);
			}
			printf("\n");
		}*/
		load_people_current_state(t);
		load_I_people_num();
		/*printf("I_num:\n");
		for (i = 0; i < NETWORK_SIZE; i++)
		{
			printf("%d ", I_people_num[i]);
		}
		printf("\n");
		printf("\n");*/
		load_people_fin_state();
	}
}


double load_I_pro_zong()//��t_test��ƽ�����
{
	int i;
	double pro_zong, zong=0;
	printf("I_pro_t:\n");
	for (i = 0; i < t_test; i++)
	{
		printf("%d:%f \n", i,I_pro_t[i]);
	}
	printf("\n");
	printf("\n");
	for (i = 0; i < t_test ; i++)
	{
		zong = I_pro_t[i] + zong;
	}
	pro_zong = zong/(double)t_test;
	return pro_zong;
}
double load_I_sus()
{
	int i;
	double a = 0, b = 0, c = 0, S = 0;
	for (i = 0; i < t_test; i++)
	{
		if (I_pro_t[i] != 0)
		{
			a = a+pow(I_pro_t[i], 2);
		}
		b = b+I_pro_t[i] ;
		c =c+ I_pro_t[i] ;
	}
	a = a / (double)t_test;
	b = b / (double)t_test;
	if (b != 0)
	{
		b = pow(b, 2);
	}
	c = c / (double)t_test;
	if (c == 0)
	{
		S = 0;
	}
	else
	{
		S = NETWORK_SIZE * (a - b) / c;
	}
	return S;
	
}
void load_write_file()
{
	int i,j;
	double pro;
	double S;
	errno_t err1;
	err1 = fopen_s(&fp4, "help2.txt", "w");
	//errno_t err2;
	//err2 = fopen_s(&fp5, "help2.txt", "w");
	if (err1 != 0)
	{
		puts("���ܴ��ļ�");
	}
	/*if (err2 != 0)
	{
		puts("���ܴ��ļ�sus");
	}*/
	for (la = -2; la <= -2; la = la + 2)//3
	{
		fprintf_s(fp4, "la=%f\n", la);
		//fprintf_s(fp5, "la=%f\n", la);
		for (pd = 1; pd <= 1; pd = pd + 0.02)
		{
			//fprintf_s(fp4, "pd=%f\n", pd);
			//fprintf_s(fp5, "pd=%f\n", pd);
			for (lambda = 0.001; lambda <= 0.001; lambda = lambda + 0.0002)//50
			{
				load_I_pro_limit();
				pro = load_I_pro_zong();
				fprintf_s(fp4, "%f ", pro);
				//S = load_I_sus();
				//fprintf_s(fp5, "%f ", S);
			}
			fprintf_s(fp4, "\n");
			//fprintf_s(fp5, "\n");
	    }
		fprintf_s(fp4, "\n");
		fprintf_s(fp4, "\n");
		fprintf_s(fp4, "\n");
		/*fprintf_s(fp5, "\n");
		fprintf_s(fp5, "\n");
		fprintf_s(fp5, "\n");*/
    }
	fclose(fp4);
	//fclose(fp5);
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