#include <stdio.h>
#include <string.h>
typedef struct node
{
	int w,s;
}node;
#define NUM 501
node sum[NUM][512];
int num[NUM][512],pre[NUM][512],cur[NUM][512];

int main()
{
	int w1,s1,d1,w2,s2,d2,w3,s3,d3,c1,c2,c3,d4,wi,si;
	int w0,s0,i,j,k,l,n,ret,curx,g1,g2,g3,t=1;
	while(scanf("%d",&n)&&n)
	{
		if(t!=1) printf("\n");
		scanf("%d%d%d",&w1,&s1,&d1);
		scanf("%d%d%d",&w2,&s2,&d2);
		scanf("%d%d%d",&w3,&s3,&d3);
		scanf("%d%d%d%d",&c1,&c2,&c3,&d4);
		sum[0][0].w=0;
		sum[0][0].s=0;
		for(i=0;i<NUM;i++)
		{
			sum[i+1][0].w=sum[i][0].w+w1;
			sum[i+1][0].s=sum[i][0].s+s1;
		}
		for(i=0;i<NUM;i++)
		for(j=0;j<NUM;j++)
		{
			sum[i][j+1].w=sum[i][j].w+w2;
			sum[i][j+1].s=sum[i][j].s+s2;
		}
		for(i=0;i<NUM;i++)
			for(j=0;j<NUM;j++)
			{pre[i][j]=-1;cur[i][j]=-1;}
		pre[0][0]=0;
		while(n--)
		{
			scanf("%d%d",&wi,&si);
			for(i=0;i<NUM;i++)
			{
				for(j=0;j<NUM;j++)
				{
					w0 = (wi-sum[i][j].w);
					s0 = (si-sum[i][j].s);
					if(w0<0 || s0<0) 
					{num[i][j]=-1;break;}
					else
					{w0=w0/w3;s0=s0/s3;num[i][j] = w0<s0?w0:s0;}
				}
				if(j==0) break;
			}
			for(i=0;i<NUM;i++)
			{
				for(j=0;j<NUM;j++)
				{
					for(k=0;k<=i;k++)
					{
						for(l=0;l<=j;l++)
						{
							if(num[k][l]>=0){
								if(pre[i-k][j-l]>=0){
								curx = pre[i-k][j-l]+num[k][l];
								if(curx>cur[i][j])
									cur[i][j]=curx;
								}
							}
							else break;
						}
						if(l==0) break;
					}
					if(k==0) break;
				}
			}
			for(i=0;i<NUM;i++)
			{
				for(j=0;j<NUM;j++)
				{
					pre[i][j]=cur[i][j];
					if(pre[i][j]<0) break;
				}
				if(j==0) break;
			}
		}
		ret = 0;
		for(i=0;i<NUM;i++)
		{
			for(j=0;j<NUM;j++)
			{
				g1 = i;g2=j;
				g3 = cur[i][j];
				if(g3>=0)
				{
					curx = 0;
					while (g1>=c1&&g2>=c2&&g3>=c3)
					{curx+=d4;g1-=c1;g2-=c2;g3-=c3;} 
					curx += g1*d1+g2*d2+g3*d3;
					if(curx > ret) ret = curx;
				}
				else break;
			}
			if(j==0) break;
		}
		printf("Case %d: %d\n",t,ret);
		t++;
	}
	return 0;
}
