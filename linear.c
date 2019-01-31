/* Monte Calo Simulation with chain energy
 for M chains(length N) in D*D*D cube */
/*�����������ۻ������������Ҫ�Ĵ�һ��ģ��ĳ�ʼ̬*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DX 32
#define DY 32
#define DZ 32
#define M 1
#define N 32
#define F0 "Parameter.dat"
#define F1 "Inistate.dat"
#define F2 "Probability.dat"
#define F3 "FreeEnergy.dat"
#define E 1 /*this is the Ec/Ep value */
#define TempNUM 0 //the chain number of tempelate
#define R 0.00 /*this is the Er/Ep value */
#define Sample_loopnum 200
#define Status_loopnum 500000
#define Total_loopnum 2000000

double A1;
double A2;
double A0;
int B;

double StaNum;
int StaCryNum[1024 + 1];

int T1,T2,T_step;
int k;
int cx, cy, cz;
int cx1, cy1, cz1;
int bjY(int x), abjY(int x, int y), tfY(int x);
int bjX(int x), abjX(int x, int y), tfX(int x);
int bjZ(int x), abjZ(int x, int y), tfZ(int x);
int fp(int a, int b, int i);
int je0(int x1, int y1, int z1, int x2, int y2, int z2), jc(int c, int i), je1(
		int x1, int y1, int z1, int x2, int y2, int z2);
int je1N(int x1, int y1, int z1, int x2, int y2, int z2,int low,int high, int i, int j);
int je0N(int x1, int y1, int z1, int x2, int y2, int z2,int low,int high, int i, int j);
int hindbend(int a, int b, int c, int x, int y, int z, int l);
double dist(int x, int y, int z)/*,calc(void)*/;
int jpv(int x, int y, int z, int s);
double ran(long t);
int clustersize();
int clustersize1();

int xx[M + 1][N + 2];
int yy[M + 1][N + 2];
int zz[M + 1][N + 2];
int xx1[M + 1][N + 2];
int yy1[M + 1][N + 2];
int zz1[M + 1][N + 2];
int ii[DX + 1][DY + 1][DZ + 1];
int jj[DX + 1][DY + 1][DZ + 1];
int ii11[DX + 1][DY + 1][DZ + 1];
int jj11[DX + 1][DY + 1][DZ + 1];
int xemax[8 + 1], xfmax[8 + 1];
int com[M + 1][N + 1];

int the[M * N + 1]; //to record bond belong to which crtstal, and 0 mean it is armorphos, and if it is negtive ,it begin to devided to two crystals
int ti[M * N + 1]; //the temp record i of the neighbor bond while find the crystal
int tj[M * N + 1]; //the temp record j of the neighbor bond while find the crystal
int tx[M * N + 1];
int ty[M * N + 1];
int tz[M * N + 1];
int crynum,crysize;
int count[3000];

char F4[60];
int len[M + 1][2 + 1];

/*****
 * The main function.
 *****/
main(void)
{
	FILE * infile; /* Input file. */
	char linebfr[50]; /* Input line buffer */
	char filename[50];
	int a,b,c,d,xc,yc,zc,xc1,yc1,zc1,i1,i2,nn=1,cycle=1,L,r,cn,n,s=0,v=0,nv,bc,bc0,bc1,a1,b1,c1,i3,j3,i4,j4,b2;
	int i,j,j1,j2,x,y,z,p,M1;
	int num=1,x1,y1,z1,cry=0,state=0,cry0=0,cry1=0,cryZ=0,x2,y2,z2,x3,y3,z3,r1,cn1; /* the run number */
	double l,l1,l2,lb=0.0,le=0.0,lc=0.0,temp,s1,ms=0.0,ms1=0.0,ms0=0.0,E1; /*the distance between neibouring beads */
	double umE[N+1];
	double FreeEnergy[N+1];
	int t;
	//time_t t;

	//initiallize the parameters:
	infile=fopen(F0,"r");
	fscanf(infile,"%i\t%i\t%i\n",&T1,&T2,&T_step);
	fclose(infile);

	//infile=fopen(F0,"a");
	//fprintf(infile,"N=%i\tA1=%.2f\tA2=%.2f\tA0=%.2f\tB=%i\t100T=%i\n",N,A1,A2,A0,B,T);
	//fclose(infile);
	//temp=T/100.0;

	//initiallize the ouput files:

	infile=fopen(F2,"w");
	fprintf(infile,"Temperature\tcycle");
	for (i=0;i<=N;i++)
	fprintf(infile,"\t%d",i);
	fprintf(infile,"\n");
	fclose(infile);

	infile=fopen(F3,"w");
	fprintf(infile,"Temperature\tcycle");
	for (i=0;i<=N;i++)
	fprintf(infile,"\t%d",i);
	fprintf(infile,"\n");
	fclose(infile);
for(t=T1;t>=T2;t-=T_step)
    {temp=t/100.0;
    cycle=1;
	/* input the initial state*/
	for (x=1;x<=DX;x++)
	for (y=1;y<=DY;y++)
	for (z=1;z<=DZ;z++)
	ii[x][y][z]=ii11[x][y][z]=0,jj[x][y][z]=jj11[x][y][z]=-1;

	StaNum=0;

	for(a=0;a<=N;a++)
	{	StaCryNum[a]=0;}
//StaCryNum[0]=1;
	for(a=0;a<=N;a++)
	{
		if(a<=B)
		{
			umE[a]=1.0*A1*(a-B)*(a-B);
			continue;
		}
		if(a>B)
		{
			umE[a]=1.0*A2*(a-B)*(a-B);
			continue;
		}
	}
	umE[0]=umE[1]+A0;

	//l=ran(-time(&t));

	infile = fopen(F1, "r" );
	if ( !infile )
	{
		printf( "Can't open input file!\n" );
		exit(1);}
	for (i=1;i<=M;i++)
	for (j=1;j<=N;j++)
	{
		fscanf(infile,"%i\t%i\t%i\n",&a,&b,&c);
		xx[i][j]=xx1[i][j]=a;
		yy[i][j]=yy1[i][j]=b;
		zz[i][j]=zz1[i][j]=c;
		a=tfX(a),b=tfY(b),c=tfZ(c);
		ii[a][b][c]=ii11[a][b][c]=i;
		jj[a][b][c]=jj11[a][b][c]=j;
		if(j==160)
		com[i][j]=2;
		else
		com[i][j]=1;
	}
	fclose(infile);

	bc0=0;
	bc0=clustersize();

	/* go through every vacancy once */
	start:
	for (k=1;k<=M*N;k++)
	{
		i=(int)(ran(1000000)*M)+1;
		if(i<=TempNUM)continue;
		j=(int)(ran(1000000)*N)+1;
		x=xx[i][j];
		y=yy[i][j];
		z=zz[i][j];
		x=tfX(x);
		y=tfY(y);
		z=tfZ(z);

		cx=x+(int)(3.0*ran(1000000))-1;
		cy=y+(int)(3.0*ran(1000000))-1;
		cz=z+(int)(3.0*ran(1000000))-1;
		cx=bjX(cx);
		cy=bjY(cy);
		cz=bjZ(cz);

		if(ii[cx][cy][cz]!=0)continue;

		/* judge the vacancy movable */
		if (cz==z&&cx!=x&&cy!=y)
		{
			if (ii[cx][y][z]==ii[x][cy][z]&&
					abs(jj[cx][y][z]-jj[x][cy][z])==1)
			continue; /*��ֹ���ԭ�еļ�*/
			goto motion;}
		if (cy==y&&cx!=x&&cz!=z)
		{
			if (ii[cx][y][z]==ii[x][y][cz]&&
					abs(jj[cx][y][z]-jj[x][y][cz])==1)
			continue;
			goto motion;}
		if (cx==x&&cy!=y&&cz!=z)
		{
			if (ii[x][cy][z]==ii[x][y][cz]&&
					abs(jj[x][cy][z]-jj[x][y][cz])==1)
			continue;
			goto motion;}
		if (cx!=x&&cy!=y&&cz!=z)
		if ((ii[x][cy][cz]==ii[cx][y][z]&&
						abs(jj[x][cy][cz]-jj[cx][y][z])==1)||
				(ii[cx][y][cz]==ii[x][cy][z]&&
						abs(jj[cx][y][cz]-jj[x][cy][z])==1)||
				(ii[cx][cy][z]==ii[x][y][cz]&&
						abs(jj[cx][cy][z]-jj[x][y][cz])==1))
		continue;

		/* motion */
		motion:

		/* the chosed bead is the head of chain*/
		if (j==1)
		{
			a=tfX(xx[i][2]),b=tfY(yy[i][2]),c=tfZ(zz[i][2]);
			l=dist(a,b,c);
			xc=abjX(xx[i][1],cx);
			yc=abjY(yy[i][1],cy);
			zc=abjZ(zz[i][1],cz);
			if (l>1.8)
			{
				p=fp(1,N,i);

				//���趯֮��
				cx1=cx;xc1=xc;
				cy1=cy;yc1=yc;
				cz1=cz;zc1=zc;

				xx1[i][0]=xc1;
				yy1[i][0]=yc1;
				zz1[i][0]=zc1;
				cx1=tfX(xx[i][p]);
				cy1=tfY(yy[i][p]);
				cz1=tfZ(zz[i][p]);
				ii11[cx1][cy1][cz1]=0;
				jj11[cx1][cy1][cz1]=-1;
				for (i1=p;i1>=1;i1--)
				{
					xx1[i][i1]=xx1[i][i1-1]; /*����һ��������ǰ�ƶ�*/
					yy1[i][i1]=yy1[i][i1-1];
					zz1[i][i1]=zz1[i][i1-1];
					jj11[tfX(xx1[i][i1])][tfY(yy1[i][i1])][tfZ(zz1[i][i1])]++;
				} /*���ڵ�i������ȡ���˵ڣ�i-1�������ε�λ�ã���ʵ�������ǵ�i�����Σ�����jjҪ+1*/
				ii11[cx][cy][cz]=i;
				jj11[cx][cy][cz]=1;
				//�˶��������

				if (p==2)
				{
					a=((xx[i][3]-xx[i][1])!=(xx[i][1]-xc)
							||(yy[i][3]-yy[i][1])!=(yy[i][1]-yc)
							||(zz[i][3]-zz[i][1])!=(zz[i][1]-zc));
					a+=((xx[i][1]-xx[i][3])!=(xx[i][3]-xx[i][4])
							||(yy[i][1]-yy[i][3])!=(yy[i][3]-yy[i][4])
							||(zz[i][1]-zz[i][3])!=(zz[i][3]-zz[i][4]));
					a--;
					a-=((xx[i][2]-xx[i][3])!=(xx[i][3]-xx[i][4])
							||(yy[i][2]-yy[i][3])!=(yy[i][3]-yy[i][4])
							||(zz[i][2]-zz[i][3])!=(zz[i][3]-zz[i][4]));
				}
				else
				{
					a=((xx[i][2]-xx[i][1])!=(xx[i][1]-xc)
							||(yy[i][2]-yy[i][1])!=(yy[i][1]-yc)
							||(zz[i][2]-zz[i][1])!=(zz[i][1]-zc));
					if (p==N)
					{
						a-=((xx[i][N]-xx[i][N-1])!=(xx[i][N-1]-xx[i][N-2])
								||(yy[i][N]-yy[i][N-1])!=(yy[i][N-1]-yy[i][N-2])
								||(zz[i][N]-zz[i][N-1])!=(zz[i][N-1]-zz[i][N-2]));
					}
					if(p==N-1)
					{
						a+=((xx[i][N]-xx[i][N-2])!=(xx[i][N-2]-xx[i][N-3])
								||(yy[i][N]-yy[i][N-2])!=(yy[i][N-2]-yy[i][N-3])
								||(zz[i][N]-zz[i][N-2])!=(zz[i][N-2]-zz[i][N-3]));
						a-=((xx[i][N-1]-xx[i][N-2])!=(xx[i][N-2]-xx[i][N-3])
								||(yy[i][N-1]-yy[i][N-2])!=(yy[i][N-2]-yy[i][N-3])
								||(zz[i][N-1]-zz[i][N-2])!=(zz[i][N-2]-zz[i][N-3]));
						a--;
					}
					if(p<N-1) a+=jc(p,i);
				}
				b=0;
				for(j1=1;j1<N;j1++)
				{
					//if(j1==N)break;
					if(com[i][j1]==2||com[i][j1+1]==2)continue;
					b+=je1N(xx1[i][j1],yy1[i][j1],zz1[i][j1],xx1[i][j1+1],yy1[i][j1+1],zz1[i][j1+1],1,p,i,j1);
					b-=je0N(xx[i][j1],yy[i][j1],zz[i][j1],xx[i][j1+1],yy[i][j1+1],zz[i][j1+1],1,p,i,j1);
				}

				bc1=0;
				bc1=clustersize1(); //��¼��֮��ᾧbond ��
				E1=umE[bc1]-umE[bc0];//�˶�ǰ��������



				l1=(E*b-a-d*R-E1)/temp;
				if (l1<0&&ran(1000000)>exp(l1))/*metropolis sampling?*/
				{
					ii11[cx][cy][cz]=0;
					jj11[cx][cy][cz]=-1;
					for(i1=p;i1>=1;i1--)
					{
						xx1[i][i1]=xx[i][i1];
						yy1[i][i1]=yy[i][i1];
						zz1[i][i1]=zz[i][i1];
						x1=tfX(xx[i][i1]);
						y1=tfY(yy[i][i1]);
						z1=tfZ(zz[i][i1]);
						ii11[x1][y1][z1]=i;
						jj11[x1][y1][z1]=i1;
					}
					continue;
				}
				xx[i][0]=xc;
				yy[i][0]=yc;
				zz[i][0]=zc;
				cx=tfX(xx[i][p]);
				cy=tfY(yy[i][p]);
				cz=tfZ(zz[i][p]);
				ii[cx][cy][cz]=0;
				jj[cx][cy][cz]=-1;
				for (i1=p;i1>=1;i1--)
				{
					xx[i][i1]=xx[i][i1-1]; /*����һ��������ǰ�ƶ�*/
					yy[i][i1]=yy[i][i1-1];
					zz[i][i1]=zz[i][i1-1];
					jj[tfX(xx[i][i1])][tfY(yy[i][i1])][tfZ(zz[i][i1])]++;
				} /*���ڵ�i������ȡ���˵ڣ�i-1�������ε�λ�ã���ʵ�������ǵ�i�����Σ�����jjҪ+1*/
				a=tfX(xx[i][1]),b=tfY(yy[i][1]),c=tfZ(zz[i][1]);
				ii[a][b][c]=i;
				jj[a][b][c]=1;
				bc0=bc1;
				continue;
			}
			if (hindbend(a,b,c,x,y,z,(int)(l*l+0.2))==1) /*hindbend��������ֵ��1�Ļ���˵���������еļ����ǲ�������˶�*/
			continue;

			cx1=cx;xc1=xc;
			cy1=cy;yc1=yc;
			cz1=cz;zc1=zc;

			ii11[cx1][cy1][cz1]=i;
			jj11[cx1][cy1][cz1]=j;
			cx1=tfX(xx1[i][j]);
			cy1=tfY(yy1[i][j]);
			cz1=tfZ(zz1[i][j]);
			ii11[cx1][cy1][cz1]=0;
			jj11[cx1][cy1][cz1]=-1;
			xx1[i][j]=xc1;
			yy1[i][j]=yc1;
			zz1[i][j]=zc1;

			a=((xx[i][3]-xx[i][2])!=(xx[i][2]-xc)
					||(yy[i][3]-yy[i][2])!=(yy[i][2]-yc)
					||(zz[i][3]-zz[i][2])!=(zz[i][2]-zc));
			a-=((xx[i][3]-xx[i][2])!=(xx[i][2]-xx[i][1])
					||(yy[i][3]-yy[i][2])!=(yy[i][2]-yy[i][1])
					||(zz[i][3]-zz[i][2])!=(zz[i][2]-zz[i][1]));

			b=0;
			b+=je1N(xx1[i][j],yy1[i][j],zz1[i][j],xx1[i][2],yy1[i][2],zz1[i][2],1,1,i,1);
			b-=je0N(xx[i][j],yy[i][j],zz[i][j],xx[i][2],yy[i][2],zz[i][2],1,1,i,1);
			b*=2;

			//���趯֮��

			//�˶��������
			bc1=0;
			bc1=clustersize1(); //��¼��֮��ᾧbond ��
			E1=umE[bc1]-umE[bc0];//�˶�ǰ��������



			l1=(E*b-a-E1)/temp;
			if (l1<0&&ran(1000000)>exp(l1))
			{
				xx1[i][j]=xx[i][j];
				yy1[i][j]=yy[i][j];
				zz1[i][j]=zz[i][j];
				ii11[cx][cy][cz]=0;
				jj11[cx][cy][cz]=-1;
				ii11[x][y][z]=i;
				jj11[x][y][z]=j;
				continue;
			}
			xx[i][1]=xc;
			yy[i][1]=yc;
			zz[i][1]=zc;
			ii[cx][cy][cz]=i;
			jj[cx][cy][cz]=1;
			cx=x;
			cy=y;
			cz=z;
			ii[x][y][z]=0;
			jj[x][y][z]=-1;
			bc0=bc1;
			continue;
		}
		/* the chosed bead is the tail of chain */
		if (j==N)
		{
			a=tfX(xx[i][N-1]),b=tfY(yy[i][N-1]),c=tfZ(zz[i][N-1]);
			l=dist(a,b,c);
			xc=abjX(xx[i][N],cx);
			yc=abjY(yy[i][N],cy);
			zc=abjZ(zz[i][N],cz);
			if (l>1.8)
			{
				p=fp(N,1,i);

				//���趯֮��
				cx1=cx;xc1=xc;
				cy1=cy;yc1=yc;
				cz1=cz;zc1=zc;

				xx1[i][N+1]=xc1;
				yy1[i][N+1]=yc1;
				zz1[i][N+1]=zc1;
				cx1=tfX(xx[i][p]);
				cy1=tfY(yy[i][p]);
				cz1=tfZ(zz[i][p]);
				ii11[cx1][cy1][cz1]=0;
				jj11[cx1][cy1][cz1]=-1;
				for (i1=p;i1<=N;i1++)
				{
					xx1[i][i1]=xx1[i][i1+1];
					yy1[i][i1]=yy1[i][i1+1];
					zz1[i][i1]=zz1[i][i1+1];
					jj11[tfX(xx1[i][i1])][tfY(yy1[i][i1])][tfZ(zz1[i][i1])]--;
				}
				ii11[cx][cy][cz]=i;
				jj11[cx][cy][cz]=N;

				if (p==N-1)
				{
					a=((xx[i][N-2]-xx[i][N])!=(xx[i][N]-xc)
							||(yy[i][N-2]-yy[i][N])!=(yy[i][N]-yc)
							||(zz[i][N-2]-zz[i][N])!=(zz[i][N]-zc));
					a+=((xx[i][N]-xx[i][N-2])!=(xx[i][N-2]-xx[i][N-3])
							||(yy[i][N]-yy[i][N-2])!=(yy[i][N-2]-yy[i][N-3])
							||(zz[i][N]-zz[i][N-2])!=(zz[i][N-2]-zz[i][N-3]));
					a--;
					a-=((xx[i][N-1]-xx[i][N-2])!=(xx[i][N-2]-xx[i][N-3])
							||(yy[i][N-1]-yy[i][N-2])!=(yy[i][N-2]-yy[i][N-3])
							||(zz[i][N-1]-zz[i][N-2])!=(zz[i][N-2]-zz[i][N-3]));}
				else
				{
					a=((xx[i][N-1]-xx[i][N])!=(xx[i][N]-xc)
							||(yy[i][N-1]-yy[i][N])!=(yy[i][N]-yc)
							||(zz[i][N-1]-zz[i][N])!=(zz[i][N]-zc));
					if (p==1)
					{
						a-=((xx[i][1]-xx[i][2])!=(xx[i][2]-xx[i][3])
								||(yy[i][1]-yy[i][2])!=(yy[i][2]-yy[i][3])
								||(zz[i][1]-zz[i][2])!=(zz[i][2]-zz[i][3]));
					}
					if(p==2)
					{
						a+=((xx[i][1]-xx[i][3])!=(xx[i][3]-xx[i][4])
								||(yy[i][1]-yy[i][3])!=(yy[i][3]-yy[i][4])
								||(zz[i][1]-zz[i][3])!=(zz[i][3]-zz[i][4]));
						a-=((xx[i][2]-xx[i][3])!=(xx[i][3]-xx[i][4])
								||(yy[i][2]-yy[i][3])!=(yy[i][3]-yy[i][4])
								||(zz[i][2]-zz[i][3])!=(zz[i][3]-zz[i][4]));
						a--;
					}
					if(p>2) a+=jc(p,i);
				}

				b=0;
				for(j1=1;j1<N;j1++)
				{
					//if(j1==0)continue;
					if(com[i][j1]==2||com[i][j1+1]==2)continue;
					b+=je1N(xx1[i][j1],yy1[i][j1],zz1[i][j1],xx1[i][j1+1],yy1[i][j1+1],zz1[i][j1+1],p-1,N,i,j1);
					b-=je0N(xx[i][j1],yy[i][j1],zz[i][j1],xx[i][j1+1],yy[i][j1+1],zz[i][j1+1],p-1,N,i,j1);
				}

				//�˶��������
				bc1=0;
				bc1=clustersize1(); //��¼��֮��ᾧbond ��
				E1=umE[bc1]-umE[bc0];//�˶�ǰ��������


				l1=(E*b-a-d*R-E1)/temp;
				if (l1<0&&ran(1000000)>exp(l1))
				{
					ii11[cx][cy][cz]=0;
					jj11[cx][cy][cz]=-1;
					for(i1=p;i1<=N;i1++)
					{
						xx1[i][i1]=xx[i][i1];
						yy1[i][i1]=yy[i][i1];
						zz1[i][i1]=zz[i][i1];
						x1=tfX(xx[i][i1]);
						y1=tfY(yy[i][i1]);
						z1=tfZ(zz[i][i1]);
						ii11[x1][y1][z1]=i;
						jj11[x1][y1][z1]=i1;
					}
					continue;
				}
				xx[i][N+1]=xc;
				yy[i][N+1]=yc;
				zz[i][N+1]=zc;
				cx=tfX(xx[i][p]);
				cy=tfY(yy[i][p]);
				cz=tfZ(zz[i][p]);
				ii[cx][cy][cz]=0;
				jj[cx][cy][cz]=-1;
				for (i1=p;i1<=N;i1++)
				{
					xx[i][i1]=xx[i][i1+1];
					yy[i][i1]=yy[i][i1+1];
					zz[i][i1]=zz[i][i1+1];
					jj[tfX(xx[i][i1])][tfY(yy[i][i1])][tfZ(zz[i][i1])]--;}
				a=tfX(xx[i][N]),b=tfY(yy[i][N]),c=tfZ(zz[i][N]);
				ii[a][b][c]=i;
				jj[a][b][c]=N;
				bc0=bc1;
				continue;
			}
			if (hindbend(a,b,c,x,y,z,(int)(l*l+0.2))==1)
			continue;

			//���趯֮��
			cx1=cx;xc1=xc;
			cy1=cy;yc1=yc;
			cz1=cz;zc1=zc;

			ii11[cx1][cy1][cz1]=i;
			jj11[cx1][cy1][cz1]=j;
			cx1=tfX(xx1[i][j]);
			cy1=tfY(yy1[i][j]);
			cz1=tfZ(zz1[i][j]);
			ii11[cx1][cy1][cz1]=0;
			jj11[cx1][cy1][cz1]=-1;
			xx1[i][j]=xc1;
			yy1[i][j]=yc1;
			zz1[i][j]=zc1;

			a=((xx[i][N-2]-xx[i][N-1])!=(xx[i][N-1]-xc)
					||(yy[i][N-2]-yy[i][N-1])!=(yy[i][N-1]-yc)
					||(zz[i][N-2]-zz[i][N-1])!=(zz[i][N-1]-zc));
			a-=((xx[i][N-2]-xx[i][N-1])!=(xx[i][N-1]-xx[i][N])
					||(yy[i][N-2]-yy[i][N-1])!=(yy[i][N-1]-yy[i][N])
					||(zz[i][N-2]-zz[i][N-1])!=(zz[i][N-1]-zz[i][N]));
			b=0;
			b+=je1N(xx1[i][j],yy1[i][j],zz1[i][j],xx1[i][N-1],yy1[i][N-1],zz1[i][N-1],N-1,N-1,i,j);
			b-=je0N(xx[i][j],yy[i][j],zz[i][j],xx[i][N-1],yy[i][N-1],zz[i][N-1],N-1,N-1,i,j);
            b*=2;

			//�˶��������
			bc1=0;
			bc1=clustersize1(); //��¼��֮��ᾧbond ��
			E1=umE[bc1]-umE[bc0];//�˶�ǰ��������



			l1=(E*b-a-E1)/temp;
			if (l1<0&&ran(1000000)>exp(l1))
			{
				xx1[i][j]=xx[i][j];
				yy1[i][j]=yy[i][j];
				zz1[i][j]=zz[i][j];
				ii11[cx][cy][cz]=0;
				jj11[cx][cy][cz]=-1;
				ii11[x][y][z]=i;
				jj11[x][y][z]=j;
				continue;
			}
			xx[i][N]=xc;
			yy[i][N]=yc;
			zz[i][N]=zc;
			ii[cx][cy][cz]=i;
			jj[cx][cy][cz]=N;
			cx=x;
			cy=y;
			cz=z;
			ii[x][y][z]=0;
			jj[x][y][z]=-1;
			bc0=bc1;
			continue;
		}

		/* the chosed bead is in middle of chain */
		a=tfX(xx[i][j-1]),b=tfY(yy[i][j-1]),c=tfZ(zz[i][j-1]);
		l=dist(a,b,c);
		xc=tfX(xx[i][j+1]),yc=tfY(yy[i][j+1]),zc=tfZ(zz[i][j+1]);
		l1=dist(xc,yc,zc);
		if (l>1.8&&l1>1.8) /*��Ҫͬʱ�϶��ڽ��������㣬������*/
		continue;
		if (l<1.8&&l1<1.8)
		{ /*���������㶼����Ҫ��*/
			if (hindbend(a,b,c,x,y,z,(int)(l*l+0.2))==1)
			continue;
			if (hindbend(xc,yc,zc,x,y,z,(int)(l1*l1+0.2))==1)
			continue;
			xc=abjX(xx[i][j],cx);
			yc=abjY(yy[i][j],cy);
			zc=abjZ(zz[i][j],cz);

			//���趯֮��
			cx1=cx;xc1=xc;
			cy1=cy;yc1=yc;
			cz1=cz;zc1=zc;

			ii11[cx1][cy1][cz1]=i;
			jj11[cx1][cy1][cz1]=j;
			cx1=tfX(xx1[i][j]);
			cy1=tfY(yy1[i][j]);
			cz1=tfZ(zz1[i][j]);
			ii11[cx1][cy1][cz1]=0;
			jj11[cx1][cy1][cz1]=-1;
			xx1[i][j]=xc1;
			yy1[i][j]=yc1;
			zz1[i][j]=zc1;

			a=0;

			if (j<N-1)
			{
				a+=((xx[i][j+2]-xx[i][j+1])!=(xx[i][j+1]-xc)
						||(yy[i][j+2]-yy[i][j+1])!=(yy[i][j+1]-yc)
						||(zz[i][j+2]-zz[i][j+1])!=(zz[i][j+1]-zc));
				a-=((xx[i][j+2]-xx[i][j+1])!=(xx[i][j+1]-xx[i][j])
						||(yy[i][j+2]-yy[i][j+1])!=(yy[i][j+1]-yy[i][j])
						||(zz[i][j+2]-zz[i][j+1])!=(zz[i][j+1]-zz[i][j]));
			}
			if (j>2)
			{
				a+=((xx[i][j-2]-xx[i][j-1])!=(xx[i][j-1]-xc)
						||(yy[i][j-2]-yy[i][j-1])!=(yy[i][j-1]-yc)
						||(zz[i][j-2]-zz[i][j-1])!=(zz[i][j-1]-zc));
				a-=((xx[i][j-2]-xx[i][j-1])!=(xx[i][j-1]-xx[i][j])
						||(yy[i][j-2]-yy[i][j-1])!=(yy[i][j-1]-yy[i][j])
						||(zz[i][j-2]-zz[i][j-1])!=(zz[i][j-1]-zz[i][j]));
			}
			a+=((xx[i][j+1]-xc)!=(xc-xx[i][j-1])
					||(yy[i][j+1]-yc)!=(yc-yy[i][j-1])
					||(zz[i][j+1]-zc)!=(zc-zz[i][j-1]));
			a-=((xx[i][j+1]-xx[i][j])!=(xx[i][j]-xx[i][j-1])
					||(yy[i][j+1]-yy[i][j])!=(yy[i][j]-yy[i][j-1])
					||(zz[i][j+1]-zz[i][j])!=(zz[i][j]-zz[i][j-1]));
			b=0;
			b+=je1N(xx1[i][j],yy1[i][j],zz1[i][j],xx1[i][j+1],yy1[i][j+1],zz1[i][j+1],j,j+1,i,j);
			b-=je0N(xx[i][j],yy[i][j],zz[i][j],xx[i][j+1],yy[i][j+1],zz[i][j+1],j,j+1,i,j);
			b+=je1N(xx1[i][j],yy1[i][j],zz1[i][j],xx1[i][j-1],yy1[i][j-1],zz1[i][j-1],j-1,j,i,j-1);
			b-=je0N(xx[i][j],yy[i][j],zz[i][j],xx[i][j-1],yy[i][j-1],zz[i][j-1],j-1,j,i,j-1);
			b*=2;

			//�˶��������
			bc1=0;
			bc1=clustersize1(); //��¼��֮��ᾧbond ��
			E1=umE[bc1]-umE[bc0];//�˶�ǰ��������



			l1=(E*b-a-E1)/temp;
			if (l1<0&&ran(1000000)>exp(l1))
			{
				xx1[i][j]=xx[i][j];
				yy1[i][j]=yy[i][j];
				zz1[i][j]=zz[i][j];
				ii11[cx][cy][cz]=0;
				jj11[cx][cy][cz]=-1;
				ii11[x][y][z]=i;
				jj11[x][y][z]=j;
				continue;
			}
			ii[cx][cy][cz]=i;
			jj[cx][cy][cz]=j;
			cx=tfX(xx[i][j]);
			cy=tfY(yy[i][j]);
			cz=tfZ(zz[i][j]);
			ii[cx][cy][cz]=0;
			jj[cx][cy][cz]=-1;
			xx[i][j]=xc;
			yy[i][j]=yc;
			zz[i][j]=zc;
			bc0=bc1;
			continue;}
		if (l>1.8&&l1<1.8)
		{ /*��ߵĵ�Ҫ���Ŷ�*/
			if (hindbend(xc,yc,zc,x,y,z,(int)(l1*l1+0.2))==1)
			continue;
			p=fp(j,1,i);
			xc=abjX(xx[i][j],cx);
			yc=abjY(yy[i][j],cy);
			zc=abjZ(zz[i][j],cz);

			//���趯֮��
			cx1=cx;xc1=xc;
			cy1=cy;yc1=yc;
			cz1=cz;zc1=zc;

			ii11[cx1][cy1][cz1]=i;
			jj11[cx1][cy1][cz1]=j;
			cx1=tfX(xx[i][p]);
			cy1=tfY(yy[i][p]);
			cz1=tfZ(zz[i][p]);
			ii11[cx1][cy1][cz1]=0;
			jj11[cx1][cy1][cz1]=-1;
			for(i1=p;i1!=j;i1=j2)
			{
				j2=i1+1;
				xx1[i][i1]=xx1[i][j2];
				yy1[i][i1]=yy1[i][j2];
				zz1[i][i1]=zz1[i][j2];
				jj11[tfX(xx1[i][i1])][tfY(yy1[i][i1])][tfZ(zz1[i][i1])]=i1;
			}
			xx1[i][j]=xc1;
			yy1[i][j]=yc1;
			zz1[i][j]=zc1;

			a=0;
			if (p==j-1)
			{
				if (j<N-1)
				{
					a+=((xx[i][j+2]-xx[i][j+1])!=(xx[i][j+1]-xc)
							||(yy[i][j+2]-yy[i][j+1])!=(yy[i][j+1]-yc)
							||(zz[i][j+2]-zz[i][j+1])!=(zz[i][j+1]-zc));
					a-=((xx[i][j+2]-xx[i][j+1])!=(xx[i][j+1]-xx[i][j])
							||(yy[i][j+2]-yy[i][j+1])!=(yy[i][j+1]-yy[i][j])
							||(zz[i][j+2]-zz[i][j+1])!=(zz[i][j+1]-zz[i][j]));
				} // end if j<N-1
				a+=((xx[i][j]-xc)!=(xc-xx[i][j+1])
						||(yy[i][j]-yc)!=(yc-yy[i][j+1])
						||(zz[i][j]-zc)!=(zc-zz[i][j+1]));
				a-=((xx[i][j-1]-xx[i][j])!=(xx[i][j]-xx[i][j+1])
						||(yy[i][j-1]-yy[i][j])!=(yy[i][j]-yy[i][j+1])
						||(zz[i][j-1]-zz[i][j])!=(zz[i][j]-zz[i][j+1]));
				if (p!=1)
				{
					a+=((xx[i][j-2]-xx[i][j])!=(xx[i][j]-xc)
							||(yy[i][j-2]-yy[i][j])!=(yy[i][j]-yc)
							||(zz[i][j-2]-zz[i][j])!=(zz[i][j]-zc));
					a--;
					if (p!=2)
					{
						a+=((xx[i][j]-xx[i][j-2])!=(xx[i][j-2]-xx[i][j-3])
								||(yy[i][j]-yy[i][j-2])!=(yy[i][j-2]-yy[i][j-3])
								||(zz[i][j]-zz[i][j-2])!=(zz[i][j-2]-zz[i][j-3]));
						a-=((xx[i][j-1]-xx[i][j-2])!=(xx[i][j-2]-xx[i][j-3])
								||(yy[i][j-1]-yy[i][j-2])!=(yy[i][j-2]-yy[i][j-3])
								||(zz[i][j-1]-zz[i][j-2])!=(zz[i][j-2]-zz[i][j-3]));
					}
				}
			}	//end if p==j-1
			else
			{
				if (j<N-1)
				{
					a+=((xx[i][j+2]-xx[i][j+1])!=(xx[i][j+1]-xc)
							||(yy[i][j+2]-yy[i][j+1])!=(yy[i][j+1]-yc)
							||(zz[i][j+2]-zz[i][j+1])!=(zz[i][j+1]-zc));
					a-=((xx[i][j+2]-xx[i][j+1])!=(xx[i][j+1]-xx[i][j])
							||(yy[i][j+2]-yy[i][j+1])!=(yy[i][j+1]-yy[i][j])
							||(zz[i][j+2]-zz[i][j+1])!=(zz[i][j+1]-zz[i][j]));
				}
				a+=((xx[i][j-1]-xx[i][j])!=(xx[i][j]-xc)
						||(yy[i][j-1]-yy[i][j])!=(yy[i][j]-yc)
						||(zz[i][j-1]-zz[i][j])!=(zz[i][j]-zc));
				a+=((xx[i][j]-xc)!=(xc-xx[i][j+1])
						||(yy[i][j]-yc)!=(yc-yy[i][j+1])
						||(zz[i][j]-zc)!=(zc-zz[i][j+1]));
				a-=((xx[i][j-1]-xx[i][j])!=(xx[i][j]-xx[i][j+1])
						||(yy[i][j-1]-yy[i][j])!=(yy[i][j]-yy[i][j+1])
						||(zz[i][j-1]-zz[i][j])!=(zz[i][j]-zz[i][j+1]));
				if (p==1)
				{
					a-=((xx[i][1]-xx[i][2])!=(xx[i][2]-xx[i][3])
							||(yy[i][1]-yy[i][2])!=(yy[i][2]-yy[i][3])
							||(zz[i][1]-zz[i][2])!=(zz[i][2]-zz[i][3]));
				}
				if(p==2)
				{
					a+=((xx[i][1]-xx[i][3])!=(xx[i][3]-xx[i][4])
							||(yy[i][1]-yy[i][3])!=(yy[i][3]-yy[i][4])
							||(zz[i][1]-zz[i][3])!=(zz[i][3]-zz[i][4]));
					a-=((xx[i][2]-xx[i][3])!=(xx[i][3]-xx[i][4])
							||(yy[i][2]-yy[i][3])!=(yy[i][3]-yy[i][4])
							||(zz[i][2]-zz[i][3])!=(zz[i][3]-zz[i][4]));
					a--;
				}
				if(p>2) a+=jc(p,i);
			}

			b=0;
			for(j1=1;j1<N;j1++)
			{
				//if(j1==0)continue;
				if(com[i][j1]==2||com[i][j1+1]==2)continue;
				b+=je1N(xx1[i][j1],yy1[i][j1],zz1[i][j1],xx1[i][j1+1],yy1[i][j1+1],zz1[i][j1+1],p-1,j,i,j1);
				b-=je0N(xx[i][j1],yy[i][j1],zz[i][j1],xx[i][j1+1],yy[i][j1+1],zz[i][j1+1],p-1,j,i,j1);
			}

			//�˶��������
			bc1=0;
		    bc1=clustersize1(); //��¼��֮��ᾧbond ��
			E1=umE[bc1]-umE[bc0];//�˶�ǰ��������



			l1=(E*b-a-d*R-E1)/temp;
			if (l1<0&&ran(1000000)>exp(l1))
			{
				ii11[cx][cy][cz]=0;
				jj11[cx][cy][cz]=-1;
				j1=j+1;
				if(j1>N)j1-=N;
				for(i1=p;i1!=j1;i1=j2)
				{
					j2=i1+1;
					if(j2>N)j2-=N;
					xx1[i][i1]=xx[i][i1];
					yy1[i][i1]=yy[i][i1];
					zz1[i][i1]=zz[i][i1];
					x1=tfX(xx[i][i1]);
					y1=tfY(yy[i][i1]);
					z1=tfZ(zz[i][i1]);
					ii11[x1][y1][z1]=i;
					jj11[x1][y1][z1]=i1;
				}
				continue;
			}
			ii[cx][cy][cz]=i; /*ԭ��Ϊ�ո�ĵط��Ѿ������ϵĵ���䣬����Ҫ������������������*/
			jj[cx][cy][cz]=j;
			cx=tfX(xx[i][p]); /*��p�����㶯�ĵ㻹��û���ĵ㣿������ǡ���㶯�ĵ㣩*/
			cy=tfY(yy[i][p]);
			cz=tfZ(zz[i][p]); /*��Ϊp�����ˣ�����ԭ��p��λ�������Ǹ��ո�*/
			ii[cx][cy][cz]=0;
			jj[cx][cy][cz]=-1;
			for (i1=p;i1<=j-1;i1++)
			{ /*��ߵ��������Ҷ�������*/
				xx[i][i1]=xx[i][i1+1];
				yy[i][i1]=yy[i][i1+1];
				zz[i][i1]=zz[i][i1+1];
				jj[tfX(xx[i][i1])][tfY(yy[i][i1])][tfZ(zz[i][i1])]--; /*����i1�������Ѿ�����ֵ���ˣ�i+1�����ˣ���ʵ�ʵ����α��ȴû�䣬��jj������ֵ��ԭ����i+1�����εı�ţ�Ҫ��һ*/
			}
			xx[i][j]=xc;
			yy[i][j]=yc;
			zz[i][j]=zc;
			bc0=bc1;
			continue;}
		if (l<1.8&&l1>1.8)
		{ /*��߲����ұ߶�*/
			if (hindbend(a,b,c,x,y,z,(int)(l*l+0.2))==1)
			continue;
			p=fp(j,N,i);
			xc=abjX(xx[i][j],cx);
			yc=abjY(yy[i][j],cy);
			zc=abjZ(zz[i][j],cz);

			//���趯֮��
			cx1=cx;xc1=xc;
			cy1=cy;yc1=yc;
			cz1=cz;zc1=zc;

			ii11[cx1][cy1][cz1]=i;
			jj11[cx1][cy1][cz1]=j;
			cx1=tfX(xx[i][p]);
			cy1=tfY(yy[i][p]);
			cz1=tfZ(zz[i][p]);
			ii11[cx1][cy1][cz1]=0;
			jj11[cx1][cy1][cz1]=-1;
			for (i1=p;i1!=j;i1=j2)
			{
				j2=i1-1;
				if(j2<1)j2+=N;
				xx1[i][i1]=xx1[i][j2];
				yy1[i][i1]=yy1[i][j2];
				zz1[i][i1]=zz1[i][j2];
				jj11[tfX(xx1[i][i1])][tfY(yy1[i][i1])][tfZ(zz1[i][i1])]=i1;
			}
			xx1[i][j]=xc1;
			yy1[i][j]=yc1;
			zz1[i][j]=zc1;

			a=0;
			if (p==j+1)
			{
				if (j>2)
				{
					a+=((xx[i][j-2]-xx[i][j-1])!=(xx[i][j-1]-xc)
							||(yy[i][j-2]-yy[i][j-1])!=(yy[i][j-1]-yc)
							||(zz[i][j-2]-zz[i][j-1])!=(zz[i][j-1]-zc));
					a-=((xx[i][j-2]-xx[i][j-1])!=(xx[i][j-1]-xx[i][j])
							||(yy[i][j-2]-yy[i][j-1])!=(yy[i][j-1]-yy[i][j])
							||(zz[i][j-2]-zz[i][j-1])!=(zz[i][j-1]-zz[i][j]));
				}
				a+=((xx[i][j]-xc)!=(xc-xx[i][j-1])
						||(yy[i][j]-yc)!=(yc-yy[i][j-1])
						||(zz[i][j]-zc)!=(zc-zz[i][j-1]));
				a-=((xx[i][j-1]-xx[i][j])!=(xx[i][j]-xx[i][j+1])
						||(yy[i][j-1]-yy[i][j])!=(yy[i][j]-yy[i][j+1])
						||(zz[i][j-1]-zz[i][j])!=(zz[i][j]-zz[i][j+1]));
				if (p!=N)
				{
					a+=((xx[i][j+2]-xx[i][j])!=(xx[i][j]-xc)
							||(yy[i][j+2]-yy[i][j])!=(yy[i][j]-yc)
							||(zz[i][j+2]-zz[i][j])!=(zz[i][j]-zc));
					a--;
					if (p!=N-1)
					{
						a+=((xx[i][j]-xx[i][j+2])!=(xx[i][j+2]-xx[i][j+3])
								||(yy[i][j]-yy[i][j+2])!=(yy[i][j+2]-yy[i][j+3])
								||(zz[i][j]-zz[i][j+2])!=(zz[i][j+2]-zz[i][j+3]));
						a-=((xx[i][j+1]-xx[i][j+2])!=(xx[i][j+2]-xx[i][j+3])
								||(yy[i][j+1]-yy[i][j+2])!=(yy[i][j+2]-yy[i][j+3])
								||(zz[i][j+1]-zz[i][j+2])!=(zz[i][j+2]-zz[i][j+3]));
					}
				}
			}
			else
			{
				if (j>2)
				{
					a+=((xx[i][j-2]-xx[i][j-1])!=(xx[i][j-1]-xc)
							||(yy[i][j-2]-yy[i][j-1])!=(yy[i][j-1]-yc)
							||(zz[i][j-2]-zz[i][j-1])!=(zz[i][j-1]-zc));
					a-=((xx[i][j-2]-xx[i][j-1])!=(xx[i][j-1]-xx[i][j])
							||(yy[i][j-2]-yy[i][j-1])!=(yy[i][j-1]-yy[i][j])
							||(zz[i][j-2]-zz[i][j-1])!=(zz[i][j-1]-zz[i][j]));
				}
				a+=((xx[i][j+1]-xx[i][j])!=(xx[i][j]-xc)
						||(yy[i][j+1]-yy[i][j])!=(yy[i][j]-yc)
						||(zz[i][j+1]-zz[i][j])!=(zz[i][j]-zc));
				a+=((xx[i][j]-xc)!=(xc-xx[i][j-1])
						||(yy[i][j]-yc)!=(yc-yy[i][j-1])
						||(zz[i][j]-zc)!=(zc-zz[i][j-1]));
				a-=((xx[i][j-1]-xx[i][j])!=(xx[i][j]-xx[i][j+1])
						||(yy[i][j-1]-yy[i][j])!=(yy[i][j]-yy[i][j+1])
						||(zz[i][j-1]-zz[i][j])!=(zz[i][j]-zz[i][j+1]));
				if (p==N)
				{
					a-=((xx[i][N]-xx[i][N-1])!=(xx[i][N-1]-xx[i][N-2])
							||(yy[i][N]-yy[i][N-1])!=(yy[i][N-1]-yy[i][N-2])
							||(zz[i][N]-zz[i][N-1])!=(zz[i][N-1]-zz[i][N-2]));
				}
				if(p==N-1)
				{
					a+=((xx[i][N]-xx[i][N-2])!=(xx[i][N-2]-xx[i][N-3])
							||(yy[i][N]-yy[i][N-2])!=(yy[i][N-2]-yy[i][N-3])
							||(zz[i][N]-zz[i][N-2])!=(zz[i][N-2]-zz[i][N-3]));
					a-=((xx[i][N-1]-xx[i][N-2])!=(xx[i][N-2]-xx[i][N-3])
							||(yy[i][N-1]-yy[i][N-2])!=(yy[i][N-2]-yy[i][N-3])
							||(zz[i][N-1]-zz[i][N-2])!=(zz[i][N-2]-zz[i][N-3]));
					a--;
				}
				if(p<N-1) a+=jc(p,i);
			}

			b=0;
			for(j1=1;j1<N;j1++)
			{
				//if(j1==N)break;
				if(com[i][j1]==2||com[i][j1+1]==2)continue;
				b+=je1N(xx1[i][j1],yy1[i][j1],zz1[i][j1],xx1[i][j1+1],yy1[i][j1+1],zz1[i][j1+1],j-1,p,i,j1);
				b-=je0N(xx[i][j1],yy[i][j1],zz[i][j1],xx[i][j1+1],yy[i][j1+1],zz[i][j1+1],j-1,p,i,j1);
			}

			//�˶��������
			bc1=0;
			bc1=clustersize1(); //��¼��֮��ᾧbond ��
			E1=umE[bc1]-umE[bc0];//�˶�ǰ��������



			l1=(E*b-a-d*R-E1)/temp;
			if (l1<0&&ran(1000000)>exp(l1))
			{
				ii11[cx][cy][cz]=0;
				jj11[cx][cy][cz]=-1;
				j1=j-1;
				if(j1<1)j1+=N;
				for(i1=p;i1!=j1;i1=j2)
				{
					j2=i1-1;
					if(j2<1)j2+=N;
					xx1[i][i1]=xx[i][i1];
					yy1[i][i1]=yy[i][i1];
					zz1[i][i1]=zz[i][i1];
					x1=tfX(xx[i][i1]);
					y1=tfY(yy[i][i1]);
					z1=tfZ(zz[i][i1]);
					ii11[x1][y1][z1]=i;
					jj11[x1][y1][z1]=i1;
				}
				continue;
			}
			ii[cx][cy][cz]=i;
			jj[cx][cy][cz]=j;
			cx=tfX(xx[i][p]);
			cy=tfY(yy[i][p]);
			cz=tfZ(zz[i][p]);
			ii[cx][cy][cz]=0;
			jj[cx][cy][cz]=-1;
			for (i1=p;i1>=j+1;i1--)
			{
				xx[i][i1]=xx[i][i1-1];
				yy[i][i1]=yy[i][i1-1];
				zz[i][i1]=zz[i][i1-1];
				jj[tfX(xx[i][i1])][tfY(yy[i][i1])][tfZ(zz[i][i1])]++;}
			xx[i][j]=xc;
			yy[i][j]=yc;
			zz[i][j]=zc;
			bc0=bc1;
			continue;}
	}

	/* display the run number and time */
	cycle++;

	if (cycle%Sample_loopnum==0)
	{
		StaNum++;

		for (i=1;i<=M;i++)
		{
			bc=0;
			bc=clustersize();

			StaCryNum[bc]++;
            /*
            if(bc==12){
			sprintf(filename,"crystalnun%istate%i.ord",bc,state++);
			if(state==100)
			state=0;
			infile=fopen(filename,"w");
			for(i=1;i<=M;i++)
			for(j=1;j<=N;j++)
			fprintf(infile,"%i\t%i\t%i\t%i\n",xx[i][j],yy[i][j],zz[i][j],com[i][j]);
			fclose(infile);}
			//check
			if(bc!=bc0)printf("wrong_bc\n");
			for(i=1;i<=M;i++)
			for(j=1;j<=N;j++)
			{
				if(xx[i][j]!=xx1[i][j]||yy[i][j]!=yy1[i][j]||zz[i][j]!=zz1[i][j])
				printf("wrong_xx\n");
			}
			for(x=1;x<=DX;x++)
			for(y=1;y<=DY;y++)
			for(z=1;z<=DZ;z++)
			{
				if(ii[x][y][z]!=ii11[x][y][z]||jj[x][y][z]!=jj11[x][y][z])
				printf("wrong_ii\n");
			}
			*/

		}

		if(cycle%Status_loopnum==0)
		{
            printf("%i\n",cycle);
			infile = fopen(F2,"a");
			fprintf(infile,"%.2f\t%i",temp,cycle);
			for(a=0;a<=N;a++)
			{
				l=1.0*StaCryNum[a]/StaNum;
				fprintf(infile,"\t%.6f",l);
				FreeEnergy[a]=-log(l);
			}
			fprintf(infile,"\n");
			fclose(infile);

			//count off the umbrella sampling
			/*
			for(a=0;a<=N;a++)
			{
				FreeEnergy[a]=FreeEnergy[a]-umE[a]/temp;
			}
			*/

			//Set the Free energy 0 point
			l=FreeEnergy[0];
			if (l!=-log(0))
			for(a=0;a<=N;a++)
			{
				FreeEnergy[a]=FreeEnergy[a]-l;
			}

			infile = fopen(F3,"a");
			fprintf(infile,"%.2f\t%i",temp,cycle);
			for(a=0;a<=N;a++)
			{
				fprintf(infile,"\t%.6f",FreeEnergy[a]);
			}
			fprintf(infile,"\n");

			/* fprintf(infile,"%.2f\t%i",temp,cycle);
			 for(a=0;a<=N;a++)
			 {
			 fprintf(infile,"\t%.6f",-log(1.0*intra[a]/intrapara[0])*temp-(a-B)*(a-B)*A);
			 }
			 fprintf(infile,"\n");  */
			fclose(infile);

			/*
			 sprintf(filename,"ring_cycle=%i_100T=%i.ord",cycle,(int)(temp*100));
			 infile=fopen(filename,"w");
			 for (i=1;i<=M;i++)
			 for (j=1;j<=N;j++)
			 fprintf(infile,"%i	%i	%i\n",xx[i][j],yy[i][j],zz[i][j]);
			 fclose(infile);
			 */
		}
	}

	if (cycle>=Total_loopnum) continue;
	goto start;
}
}
/*	if (num==16)
 {
 if (cycle%200==0)
 {
 s1=0;//����tʱ�������t=0ʱ�̶�Ӧ�ĸ�������Ԫ���� �仯�ľ���ֵ
 for(i=1;i<=M;i++)
 {
 xc=0;yc=0;zc=0;
 for(j=1;j<=N;j++)
 {
 xc=(xx[i][j]-xx0[i][j])*(xx[i][j]-xx0[i][j]);
 yc=(yy[i][j]-yy0[i][j])*(yy[i][j]-yy0[i][j]);
 zc=(zz[i][j]-zz0[i][j])*(zz[i][j]-zz0[i][j]);
 s1+=1.0*(xc+yc+zc);
 }

 }

 infile = fopen(F2,"a");
 fprintf(infile,"%i\t%.5f\n",cycle,s1/M/N);
 fclose(infile);


 sprintf(filename,"scale_chain_cycle=%i_100T=%i.ord",cycle,(int)(temp*100));
 infile=fopen(filename,"w");
 for (i=1;i<=M;i++)
 for (j=1;j<=N;j++)
 fprintf(infile,"%i	%i	%i\n",xx[i][j],yy[i][j],zz[i][j]);
 fclose(infile);


 }
 cycle++;
 num=0;
 }
 if (num>=20000000) exit(0);
 num++;
 goto start;
 }
 */

/* change the position to the initial cell */
int tfY(int x)
{
	if (x < 0)
		return (x + (1 - x / DY) * DY); /*??????????*/
	if (x % DY == 0)
		return (DY);
	return (x % DY);
}

int tfZ(int x)
{
	if (x < 0)
		return (x + (1 - x / DZ) * DZ); /*??????????*/
	if (x % DZ == 0)
		return (DZ);
	return (x % DZ);
}

int tfX(int x)
{
	if (x < 0)
		return (x + (1 - x / DX) * DX);
	if (x % DX == 0)
		return (DX);
	return (x % DX);
}

/* judge the range condition */
int bjY(int x)
{
	if (x > DY)
		x -= DY;
	else if (x < 1) /*XС����Ļ�Ҳ���Է���Ԫ��*/
		x += DY;
	return (x);
}

int bjZ(int x)
{
	if (x > DZ)
		x -= DZ;
	else if (x < 1) /*XС����Ļ�Ҳ���Է���Ԫ��*/
		x += DZ;
	return (x);
}

int bjX(int x)
{
	if (x > DX)
		x -= DX;
	else if (x < 1)
		x += DX;
	return (x);
}

/* judge the cell position */
int abjY(int x, int y)
{
	if (x < 1)
	{
		if (x % DY == 0 && y == 1)
			return (y + x / DY * DY);
		if (x % DY == 1 - DY && y == DY)
			return (y + (x / DY - 2) * DY);
		return (y + (x / DY - 1) * DY);
	}
	if (x % DY == 1 && y == DY)
		return (y + (x / DY - 1) * DY);
	if (x % DY == 0)
		if (y == DY || y == DY - 1)
			return (y + (x - 1) / DY * DY);
	return (y + x / DY * DY);
}

int abjZ(int x, int y)
{
	if (x < 1)
	{
		if (x % DZ == 0 && y == 1)
			return (y + x / DZ * DZ);
		if (x % DZ == 1 - DZ && y == DZ)
			return (y + (x / DZ - 2) * DZ);
		return (y + (x / DZ - 1) * DZ);
	}
	if (x % DZ == 1 && y == DZ)
		return (y + (x / DZ - 1) * DZ);
	if (x % DZ == 0)
		if (y == DZ || y == DZ - 1)
			return (y + (x - 1) / DZ * DZ);
	return (y + x / DZ * DZ);
}

int abjX(int x, int y)
{
	if (x < 1)
	{
		if (x % DX == 0 && y == 1)
			return (y + x / DX * DX);
		if (x % DX == 1 - DX && y == DX)
			return (y + (x / DX - 2) * DX);
		return (y + (x / DX - 1) * DX);
	}
	if (x % DX == 1 && y == DX)
		return (y + (x / DX - 1) * DX);
	if (x % DX == 0)
		if (y == DX || y == DX - 1)
			return (y + (x - 1) / DX * DX);
	return (y + x / DX * DX);
}

/* cal. the distance between vacancy and next bead */
double dist(int x, int y, int z)
{
	int xs, ys, zs;
	double xi, yi, zi;
	xs = abs(x - cx);
	if (xs > 2)
		xs = DX - xs; /*why?*/
	ys = abs(y - cy);
	if (ys > 2)
		ys = DY - ys;
	zs = abs(z - cz);
	if (zs > 2)
		zs = DZ - zs;
	xi = pow(xs, 2);
	yi = pow(ys, 2);
	zi = pow(zs, 2);
	return (sqrt(xi + yi + zi));
}

/* find the point to release */
int fp(int a, int b, int i)
{
	int i1, x, y, z, r, x1, x2, x3, y1, y2, y3, z1, z2, z3;
	if (a < b)
	{
		for (i1 = a; i1 <= b - 2; i1++)
		{
			x = xx[i][i1] - xx[i][i1 + 2];
			y = yy[i][i1] - yy[i][i1 + 2];
			z = zz[i][i1] - zz[i][i1 + 2];
			r = x * x + y * y + z * z;
			if (r > 3)
				continue;
			if (r == 1)
				return (i1 + 1);
			x1 = tfX(xx[i][i1]), x2 = tfX(xx[i][i1 + 1]), x3 = tfX(
					xx[i][i1 + 2]);
			y1 = tfY(yy[i][i1]), y2 = tfY(yy[i][i1 + 1]), y3 = tfY(
					yy[i][i1 + 2]);
			z1 = tfZ(zz[i][i1]), z2 = tfZ(zz[i][i1 + 1]), z3 = tfZ(
					zz[i][i1 + 2]);
			if (r == 2)
			{
				if (x == 0 && x1 != x2)
				{
					if (ii[x1][y3][z1] == ii[x1][y1][z3] && /*��ֹ�γɵ��¼����ԭ�����ڵľɼ�*/
					abs(jj[x1][y3][z1] - jj[x1][y1][z3]) == 1)
						continue;
					if (y2 == y1 && z2 == z3 && ii[x1][y1][z3] == ii[x2][y3][z1]
							&& abs(jj[x1][y1][z3] - jj[x2][y3][z1]) == 1) /* ��Խ��߻�����Խ��ߣ�������(ԭ����һ���ɼ���һ����ڶ�����һ���棬������Խ��ߣ����ǲ�������˶���*/
						continue;
					if (y2 == y3 && z2 == z1 && ii[x1][y3][z1] == ii[x2][y1][z3]
							&& abs(jj[x1][y3][z1] - jj[x2][y1][z3]) == 1)
						continue;
					return (i1 + 1);
				}
				if (y == 0 && y1 != y2)
				{
					if (ii[x1][y1][z3] == ii[x3][y1][z1]
							&& abs(jj[x1][y1][z3] - jj[x3][y1][z1]) == 1)
						continue;
					if (x2 == x1 && z2 == z3 && ii[x1][y1][z3] == ii[x3][y2][z1]
							&& abs(jj[x1][y1][z3] - jj[x3][y2][z1]) == 1)
						continue;
					if (x2 == x3 && z2 == z1 && ii[x3][y1][z1] == ii[x1][y2][z3]
							&& abs(jj[x3][y1][z1] - jj[x1][y2][z3]) == 1)
						continue;
					return (i1 + 1);
				}
				if (z == 0 && z1 != z2)
				{
					if (ii[x1][y3][z1] == ii[x3][y1][z1]
							&& abs(jj[x1][y3][z1] - jj[x3][y1][z1]) == 1)
						continue;
					if (x2 == x1 && y2 == y3 && ii[x1][y3][z1] == ii[x3][y1][z2]
							&& abs(jj[x1][y3][z1] - jj[x3][y1][z2]) == 1)
						continue;
					if (x2 == x3 && y2 == y1 && ii[x3][y1][z1] == ii[x1][y3][z2]
							&& abs(jj[x3][y1][z1] - jj[x1][y3][z2]) == 1)
						continue;
					return (i1 + 1);
				}
				return (i1 + 1);
			}
			if (r == 3)
			{
				if ((ii[x1][y3][z3] == ii[x3][y1][z1]
						&& abs(jj[x1][y3][z3] - jj[x3][y1][z1]) == 1)
						|| (ii[x3][y1][z3] == ii[x1][y3][z1]
								&& abs(jj[x3][y1][z3] - jj[x1][y3][z1]) == 1)
						|| (ii[x3][y3][z1] == ii[x1][y1][z3]
								&& abs(jj[x3][y3][z1] - jj[x1][y1][z3]) == 1))
					continue; /*��ֹ��Խ��߻�����*/
				return (i1 + 1);
			}
		}
		return (b);
	}
	for (i1 = a; i1 >= b + 2; i1--)
	{
		x = xx[i][i1] - xx[i][i1 - 2];
		y = yy[i][i1] - yy[i][i1 - 2];
		z = zz[i][i1] - zz[i][i1 - 2];
		r = x * x + y * y + z * z;
		if (r > 3)
			continue;
		if (r == 1)
			return (i1 - 1);
		x1 = tfX(xx[i][i1]), x2 = tfX(xx[i][i1 - 1]), x3 = tfX(xx[i][i1 - 2]);
		y1 = tfY(yy[i][i1]), y2 = tfY(yy[i][i1 - 1]), y3 = tfY(yy[i][i1 - 2]);
		z1 = tfZ(zz[i][i1]), z2 = tfZ(zz[i][i1 - 1]), z3 = tfZ(zz[i][i1 - 2]);
		if (r == 2)
		{
			if (x == 0 && x1 != x2)
			{
				if (ii[x1][y3][z1] == ii[x1][y1][z3]
						&& abs(jj[x1][y3][z1] - jj[x1][y1][z3]) == 1)
					continue;
				if (y2 == y1 && z2 == z3 && ii[x1][y1][z3] == ii[x2][y3][z1]
						&& abs(jj[x1][y1][z3] - jj[x2][y3][z1]) == 1)
					continue;
				if (y2 == y3 && z2 == z1 && ii[x1][y3][z1] == ii[x2][y1][z3]
						&& abs(jj[x1][y3][z1] - jj[x2][y1][z3]) == 1)
					continue;
				return (i1 - 1);
			}
			if (y == 0 && y1 != y2)
			{
				if (ii[x1][y1][z3] == ii[x3][y1][z1]
						&& abs(jj[x1][y1][z3] - jj[x3][y1][z1]) == 1)
					continue;
				if (x2 == x1 && z2 == z3 && ii[x1][y1][z3] == ii[x3][y2][z1]
						&& abs(jj[x1][y1][z3] - jj[x3][y2][z1]) == 1)
					continue;
				if (x2 == x3 && z2 == z1 && ii[x3][y1][z1] == ii[x1][y2][z3]
						&& abs(jj[x3][y1][z1] - jj[x1][y2][z3]) == 1)
					continue;
				return (i1 - 1);
			}
			if (z == 0 && z1 != z2)
			{
				if (ii[x1][y3][z1] == ii[x3][y1][z1]
						&& abs(jj[x1][y3][z1] - jj[x3][y1][z1]) == 1)
					continue;
				if (x2 == x1 && y2 == y3 && ii[x1][y3][z1] == ii[x3][y1][z2]
						&& abs(jj[x1][y3][z1] - jj[x3][y1][z2]) == 1)
					continue;
				if (x2 == x3 && y2 == y1 && ii[x3][y1][z1] == ii[x1][y3][z2]
						&& abs(jj[x3][y1][z1] - jj[x1][y3][z2]) == 1)
					continue;
				return (i1 - 1);
			}
			return (i1 - 1);
		}
		if (r == 3)
		{
			if ((ii[x1][y3][z3] == ii[x3][y1][z1]
					&& abs(jj[x1][y3][z3] - jj[x3][y1][z1]) == 1)
					|| (ii[x3][y1][z3] == ii[x1][y3][z1]
							&& abs(jj[x3][y1][z3] - jj[x1][y3][z1]) == 1)
					|| (ii[x3][y3][z1] == ii[x1][y1][z3]
							&& abs(jj[x3][y3][z1] - jj[x1][y1][z3]) == 1))
				continue;
			return (i1 - 1);
		}
	};
	return (b);
}

int hindbend(int a, int b, int c, int x, int y, int z, int l)
{
	long int i1, j2, p;
	j2 = ii[x][y][z] * 1000 + jj[x][y][z];
	if (l == 3)
	{
		i1 = ii[a][cy][cz] * 1000 + jj[a][cy][cz];
		p = ii[cx][b][c] * 1000 + jj[cx][b][c];
		if (i1 != j2 && p != j2 && abs(i1 - p) == 1)
			return (1); /*��������1ֵ˵���м��ᱻ���*/
		i1 = ii[cx][b][cz] * 1000 + jj[cx][b][cz];
		p = ii[a][cy][c] * 1000 + jj[a][cy][c];
		if (i1 != j2 && p != j2 && abs(i1 - p) == 1)
			return (1);
		i1 = ii[cx][cy][c] * 1000 + jj[cx][cy][c];
		p = ii[a][b][cz] * 1000 + jj[a][b][cz];
		if (i1 != j2 && p != j2 && abs(i1 - p) == 1)
			return (1);
		return (0);
	}
	if (l == 2)
	{
		if (cx == a && cx != x)
		{
			if (abs(ii[a][b][cz] - ii[a][cy][c]) == 0
					&& abs(jj[a][b][cz] - jj[a][cy][c]) == 1)
				return (1);
			if (y == b && z == cz && abs(ii[a][b][cz] - ii[x][cy][c]) == 0
					&& abs(jj[a][b][cz] - jj[x][cy][c]) == 1)
				return (1);
			if (y == cy && z == c && abs(ii[a][cy][c] - ii[x][b][cz]) == 0
					&& abs(jj[a][cy][c] - jj[x][b][cz]) == 1) /*ͬ�����ɻ���Խ��ߵĽ�����Ӱ����Խ����𣿣��𰸼��ϣ�*/
				return (1);
		}
		if (cy == b && cy != y)
		{
			if (abs(ii[a][b][cz] - ii[cx][b][c]) == 0
					&& abs(jj[a][b][cz] - jj[cx][b][c]) == 1)
				return (1);
			if (x == a && z == cz && abs(ii[a][b][cz] - ii[cx][y][c]) == 0
					&& abs(jj[a][b][cz] - jj[cx][y][c]) == 1)
				return (1);
			if (x == cx && z == c && abs(ii[cx][b][c] - ii[a][y][cz]) == 0
					&& abs(jj[cx][b][c] - jj[a][y][cz]) == 1)
				return (1);
		}
		if (cz == c && cz != z)
		{
			if (abs(ii[a][cy][c] - ii[cx][b][c]) == 0
					&& abs(jj[a][cy][c] - jj[cx][b][c]) == 1)
				return (1);
			if (x == a && y == cy && abs(ii[a][cy][c] - ii[cx][b][z]) == 0
					&& abs(jj[a][cy][c] - jj[cx][b][z]) == 1)
				return (1);
			if (x == cx && y == b && abs(ii[cx][b][c] - ii[a][cy][z]) == 0
					&& abs(jj[cx][b][c] - jj[a][cy][z]) == 1)
				return (1);
		}
		return (0);
	}
	return (0);
}

/* judge the interaction energy (if near beads neibering) *//*���ж���һ�������ڵĽ����ɣ����ǲ��������8��?)*/
int je0(int x1, int y1, int z1, int x2, int y2, int z2)
{
	int ii1, ii2, jj1, jj2;
	int xs, ys, zs, xe, ye, ze;
	int u, v, w, c = 0;
	x1 = tfX(x1);
	y1 = tfY(y1);
	z1 = tfZ(z1);
	x2 = tfX(x2);
	y2 = tfY(y2);
	z2 = tfZ(z2);
	for (u = -1; u <= 1; u++)
		for (v = -1; v <= 1; v++)
			for (w = -1; w <= 1; w++)
			{ /*����forѭ�������ɼ��ڿռ��н�����ƽ��*/
			/*	if (x1 == DX && x2 == DX && y1 == y2 && u > 0)
					continue;
				if (x1 == 1 && x2 == 1 && y1 == y2 && u < 0)
					continue;*/
				if (u == 0 && v == 0 && w == 0)
					continue;
				xs = bjX(x1 + u), ys = bjY(y1 + v), zs = bjY(z1 + w);
				xe = bjX(x2 + u), ye = bjZ(y2 + v), ze = bjZ(z2 + w);
				if (xs == x2 && ys == y2 && zs == z2)
					continue;
				if (xe == x1 && ye == y1 && ze == z1)
					continue; /*��ֹ�ҵ��ļ���ԭ���ļ�����һ��?*/
				ii1 = ii[xs][ys][zs];
				ii2 = ii[xe][ye][ze];
				jj1 = jj[xs][ys][zs];
				jj2 = jj[xe][ye][ze];
				if(com[ii1][jj1]==2||com[ii2][jj2]==2)continue;
				if (ii1 == ii2 && abs(jj1 - jj2) == 1)
					c++;
			}
	return (c);
}

int je0N(int x1, int y1, int z1, int x2, int y2, int z2, int low, int high,
		int i, int j)
{
	int ii1, ii2, jj1, jj2;
	int xs, ys, zs, xe, ye, ze;
	int u, v, w, c = 0;
	int jmin;
	x1 = tfX(x1);
	y1 = tfY(y1);
	z1 = tfZ(z1);
	x2 = tfX(x2);
	y2 = tfY(y2);
	z2 = tfZ(z2);
	for (u = -1; u <= 1; u++)
		for (v = -1; v <= 1; v++)
			for (w = -1; w <= 1; w++)
			{ /*����forѭ�������ɼ��ڿռ��н�����ƽ��*/
			/*	if (x1 == DX && x2 == DX && y1 == y2 && u > 0)
					continue;
				if (x1 == 1 && x2 == 1 && y1 == y2 && u < 0)
					continue;*/
				if (u == 0 && v == 0 && w == 0)
					continue;
				xs = bjX(x1 + u), ys = bjY(y1 + v), zs = bjY(z1 + w);
				xe = bjX(x2 + u), ye = bjZ(y2 + v), ze = bjZ(z2 + w);
				if (xs == x2 && ys == y2 && zs == z2)
					continue;
				if (xe == x1 && ye == y1 && ze == z1)
					continue; /*��ֹ�ҵ��ļ���ԭ���ļ�����һ��?*/
				ii1 = ii[xs][ys][zs];
				ii2 = ii[xe][ye][ze];
				jj1 = jj[xs][ys][zs];
				jj2 = jj[xe][ye][ze];
				if(com[ii1][jj1]==2||com[ii2][jj2]==2)continue;
				if (ii1 == ii2 && abs(jj1 - jj2) == 1)
				{
					/*
					jmin=jj1;
					if(jj1>jj2)jmin=jj2;
					if (i == ii1 && jmin >= low && jmin <= high)
					{
						if (jmin <= j)
							continue;
					}
					*/
					c++;
				}
			}
	return (c);
}
//
int je1N(int x1, int y1, int z1, int x2, int y2, int z2,int low,int high, int i, int j)
{
	int ii1, ii2, jj1, jj2;
	int xs, ys, zs, xe, ye, ze;
	int u, v, w, c = 0;
	int jmin;
	x1 = tfX(x1);
	y1 = tfY(y1);
	z1 = tfZ(z1);
	x2 = tfX(x2);
	y2 = tfY(y2);
	z2 = tfZ(z2);
	for (u = -1; u <= 1; u++)
		for (v = -1; v <= 1; v++)
			for (w = -1; w <= 1; w++)
			{ /*����forѭ�������ɼ��ڿռ��н�����ƽ��*/
			/*	if (x1 == DX && x2 == DX && y1 == y2 && u > 0)
					continue;
				if (x1 == 1 && x2 == 1 && y1 == y2 && u < 0)
					continue;*/
				if (u == 0 && v == 0 && w == 0)
					continue;
				xs = bjX(x1 + u), ys = bjY(y1 + v), zs = bjY(z1 + w);
				xe = bjX(x2 + u), ye = bjZ(y2 + v), ze = bjZ(z2 + w);
				if (xs == x2 && ys == y2 && zs == z2)
					continue;
				if (xe == x1 && ye == y1 && ze == z1)
					continue; /*��ֹ�ҵ��ļ���ԭ���ļ�����һ��?*/
				ii1 = ii11[xs][ys][zs];
				ii2 = ii11[xe][ye][ze];
				jj1 = jj11[xs][ys][zs];
				jj2 = jj11[xe][ye][ze];
				if(com[ii1][jj1]==2||com[ii2][jj2]==2)continue;
				if (ii1 == ii2 && abs(jj1 - jj2) == 1)
				{
					/*
					jmin=jj1;
					if(jj1>jj2)jmin=jj2;
					if (i == ii1 && jmin >= low && jmin <= high)
					{
						if (jmin <= j)
							continue;
					}
					*/
					c++;
				}
			}
	return (c);
}

int je1(int x1, int y1, int z1, int x2, int y2, int z2)
{
	int ii1, ii2, jj1, jj2;
	int xs, ys, zs, xe, ye, ze;
	int u, v, w, c = 0;
	x1 = tfX(x1);
	y1 = tfY(y1);
	z1 = tfZ(z1);
	x2 = tfX(x2);
	y2 = tfY(y2);
	z2 = tfZ(z2);
	for (u = -1; u <= 1; u++)
		for (v = -1; v <= 1; v++)
			for (w = -1; w <= 1; w++)
			{ /*����forѭ�������ɼ��ڿռ��н�����ƽ��*/
				/*if (x1 == DX && x2 == DX && y1 == y2 && u > 0)
					continue;
				if (x1 == 1 && x2 == 1 && y1 == y2 && u < 0)
					continue;*/
				if (u == 0 && v == 0 && w == 0)
					continue;
				xs = bjX(x1 + u), ys = bjY(y1 + v), zs = bjY(z1 + w);
				xe = bjX(x2 + u), ye = bjZ(y2 + v), ze = bjZ(z2 + w);
				if (xs == x2 && ys == y2 && zs == z2)
					continue;
				if (xe == x1 && ye == y1 && ze == z1)
					continue; /*��ֹ�ҵ��ļ���ԭ���ļ�����һ��?*/
				ii1 = ii11[xs][ys][zs];
				ii2 = ii11[xe][ye][ze];
				jj1 = jj11[xs][ys][zs];
				jj2 = jj11[xe][ye][ze];
				if(com[ii1][jj1]==2||com[ii2][jj2]==2)continue;
				if (ii1 == ii2 && abs(jj1 - jj2) == 1)
					c++;
			}
	return (c);
}
/* judge the part of conformation*/
int jc(int c, int i)
{
	int b;
	b = ((xx[i][c - 2] - xx[i][c - 1]) != (xx[i][c - 1] - xx[i][c + 1])
			|| (yy[i][c - 2] - yy[i][c - 1]) != (yy[i][c - 1] - yy[i][c + 1])
			|| (zz[i][c - 2] - zz[i][c - 1]) != (zz[i][c - 1] - zz[i][c + 1]));
	b += ((xx[i][c - 1] - xx[i][c + 1]) != (xx[i][c + 1] - xx[i][c + 2])
			|| (yy[i][c - 1] - yy[i][c + 1]) != (yy[i][c + 1] - yy[i][c + 2])
			|| (zz[i][c - 1] - zz[i][c + 1]) != (zz[i][c + 1] - zz[i][c + 2]));
	b -= ((xx[i][c - 2] - xx[i][c - 1]) != (xx[i][c - 1] - xx[i][c])
			|| (yy[i][c - 2] - yy[i][c - 1]) != (yy[i][c - 1] - yy[i][c])
			|| (zz[i][c - 2] - zz[i][c - 1]) != (zz[i][c - 1] - zz[i][c]));
	b -= ((xx[i][c + 2] - xx[i][c + 1]) != (xx[i][c + 1] - xx[i][c])
			|| (yy[i][c + 2] - yy[i][c + 1]) != (yy[i][c + 1] - yy[i][c])
			|| (zz[i][c + 2] - zz[i][c + 1]) != (zz[i][c + 1] - zz[i][c]));
	b--; /*���������jc�����ĵ��ÿɲ��룬c�������kink�㣬����(c-1),c,(c+1)�Ĺ��������Ϲ���*/
	return (b);
}

double ran(long t)
{ /*generate 2*10**18 sequential random numbers (from 0 to 1.0-1.2*exp(-7.0)) by negative seed t to initiallize,*/
	/*p282 Numerical Recipes in C 2nd ed. 1992 Cambridge University Press,London, Press WH, Teukolsky SA,Vetterling WT, Flannery BP*/
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy = 0;
	static long iv[32];
	double temp;

	if (t <= 0)
	{ /*initiallize*/
		if (0 - t < 1)
			t = 1;
		else
			t = 0 - t;
		idum2 = t;
		for (j = 39; j >= 0; j--)
		{
			k = t / 53668;
			t = 40014 * (t - k * 53668) - k * 12211;
			if (t < 0)
				t += 2147483563;
			if (j < 32)
				iv[j] = t;
		}
		iy = iv[0];
	}
	k = t / 53668;
	t = 40014 * (t - k * 53668) - k * 12211;
	if (t < 0)
		t += 2147483563;
	k = idum2 / 52774;
	idum2 = 40692 * (idum2 - k * 52774) - k * 3791;
	if (idum2 < 0)
		idum2 += 2147483399;
	j = iy / (1 + 2147483562 / 32);
	iy = iv[j] - idum2;
	iv[j] = t;
	if (iy < 0)
		iy += 2147483562;
	temp = 1.0 * iy / 2147483563;
	return (temp);
}
int clustersize()
{
	int n, i, j, i1, j1, i2, j2, x1, y1, z1, x2, y2, z2, u, v, w, xs, ys, zs,
			xe, ye, ze;
	int a;
	int bigcrysize = 0;

	for (n = 1; n <= M * N; n++)
	{
		the[n] = ti[n] = tj[n] = 0;
	}
	crynum = 0;

	for (i = 1; i <= M; i++)
		for (j = 1; j < N; j++)
		{
			/////////////////////////////////////////////////////
			if(com[i][j]==2)continue;
			////////////////////////////////////////////////////////
			j1 = j + 1;

			/////////////////////////////////////////////////////////
			if(com[i][j1]==2)continue;
			//////////////////////////////////////////////////////////
			x1 = xx[i][j];
			y1 = yy[i][j];
			z1 = zz[i][j];
			x2 = xx[i][j1];
			y2 = yy[i][j1];
			z2 = zz[i][j1];
			if (je0(x1, y1, z1, x2, y2, z2) <=5)
				continue; //if this bond is armorphos, no need
			if (the[(i - 1) * N + j] != 0)
				continue; // if this bond had been scaned, no need

			for (n = 1; n <= M * N; n++)
			{
				ti[n] = tj[n] = 0;
			}

			crysize = 1;
			crynum++;
			ti[crysize] = i;
			tj[crysize] = j;
			the[(i - 1) * N + j] = crynum;

			for (n = 1; ti[n] != 0; n++)
			{
				i1 = ti[n];
				j1 = tj[n];
				i2 = ti[n];
				j2 = tj[n] + 1;
				x1 = xx[i1][j1];
				y1 = yy[i1][j1];
				z1 = zz[i1][j1];
				x2 = xx[i2][j2];
				y2 = yy[i2][j2];
				z2 = zz[i2][j2];
				for (u = -1; u <= 1; u++)
					for (v = -1; v <= 1; v++)
						for (w = -1; w <= 1; w++)
						{
							if (u == 0 && v == 0 && w == 0)
								continue;
							xs = tfX(x1 + u), ys = tfY(y1 + v), zs = tfZ(z1 + w);
							xe = tfX(x2 + u), ye = tfY(y2 + v), ze = tfZ(z2 + w);
							if (xs == x2 && ys == y2 && zs == z2)
								continue;
							if (xe == x1 && ye == y1 && ze == z1)
								continue;
							i1 = ii[xs][ys][zs];
							i2 = ii[xe][ye][ze];
							j1 = jj[xs][ys][zs];
							j2 = jj[xe][ye][ze];
							if (i1 != i2)
								continue;
							a = abs(j1 - j2);
							if (a != 1)
								continue;
							if(com[i1][j1]==2||com[i2][j2]==2)
                                continue;
							if (je0(xs, ys, zs, xe, ye, ze) > 5)
							{

								if (j1 > j2)
									j1 = j2;
								if (the[(i1 - 1) * N + j1] != 0)
									continue;
								crysize++;
								ti[crysize] = i1;
								tj[crysize] = j1;
								the[(i1 - 1) * N + j1] = crynum;
							}

						}
			}
			count[crynum] = crysize;
			if (crysize > bigcrysize)
				bigcrysize = crysize;
		}
	return (bigcrysize);
}
int clustersize1()
{
	int n, i, j, i1, j1, i2, j2, x1, y1, z1, x2, y2, z2, u, v, w, xs, ys, zs,
			xe, ye, ze;
			int a;
	int bigcrysize = 0;

	for (n = 1; n <= M * N; n++)
	{
		the[n] = ti[n] = tj[n] = 0;
	}
	crynum = 0;

	for (i = 1; i <= M; i++)
		for (j = 1; j < N; j++)
		{
			///////////////////////////////////////////////
			if(com[i][j]==2)continue;
			////////////////////////////////////////////
			j1 = j + 1;

			///////////////////////////////////////////
			if(com[i][j1]==2)continue;
			///////////////////////////////////////


			x1 = xx1[i][j];
			y1 = yy1[i][j];
			z1 = zz1[i][j];
			x2 = xx1[i][j1];
			y2 = yy1[i][j1];
			z2 = zz1[i][j1];
			if (je1(x1, y1, z1, x2, y2, z2) <=5)
				continue; //if this bond is armorphos, no need
			if (the[(i - 1) * N + j] != 0)
				continue; // if this bond had been scaned, no need

			for (n = 1; n <= M * N; n++)
			{
				ti[n] = tj[n] = 0;
			}

			crysize = 1;
			crynum++;
			ti[crysize] = i;
			tj[crysize] = j;
			the[(i - 1) * N + j] = crynum;

			for (n = 1; ti[n] != 0; n++)
			{
				i1 = ti[n];
				j1 = tj[n];
				i2 = ti[n];
				j2 = tj[n] + 1;
				x1 = xx1[i1][j1];
				y1 = yy1[i1][j1];
				z1 = zz1[i1][j1];
				x2 = xx1[i2][j2];
				y2 = yy1[i2][j2];
				z2 = zz1[i2][j2];
				for (u = -1; u <= 1; u++)
					for (v = -1; v <= 1; v++)
						for (w = -1; w <= 1; w++)
						{
							if (u == 0 && v == 0 && w == 0)
								continue;
							xs = tfX(x1 + u), ys = tfY(y1 + v), zs = tfZ(z1 + w);
							xe = tfX(x2 + u), ye = tfY(y2 + v), ze = tfZ(z2 + w);
							if (xs == x2 && ys == y2 && zs == z2)
								continue;
							if (xe == x1 && ye == y1 && ze == z1)
								continue;
							i1 = ii11[xs][ys][zs];
							i2 = ii11[xe][ye][ze];
							j1 = jj11[xs][ys][zs];
							j2 = jj11[xe][ye][ze];
							if (i1 != i2)
								continue;
							a = abs(j1 - j2);
							if (a != 1)
								continue;
							if(com[i1][j1]==2||com[i2][j2]==2)
                                continue;
							if (je1(xs, ys, zs, xe, ye, ze) > 5)
							{

								if (j1 > j2)
									j1 = j2;
								if (the[(i1 - 1) * N + j1] != 0)
									continue;
								crysize++;
								ti[crysize] = i1;
								tj[crysize] = j1;
								the[(i1 - 1) * N + j1] = crynum;
							}

						}
			}
			count[crynum] = crysize;
			if (crysize > bigcrysize)
				bigcrysize = crysize;
		}
	return (bigcrysize);
}
