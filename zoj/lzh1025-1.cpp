#include <stdio.h>

typedef struct node
{
	int a,b;
}node;

typedef struct list
{
	struct list *next;
	int b;
}list;

node dat[5000];
list lis[5000],*ppre,*pcur;

static int cmpnode(const void *p1,const void *p2)
{
	if(((node *)p1)->a==((node *)p2)->a)
        	return ((node *)p1)->b - ((node *)p2)->b;
	else
		return ((node *)p1)->a - ((node *)p2)->a;
}
 

int main()
{
	int num,i,idx,b,count;
	scanf("%d",&count);
	while(count--)
	{
		scanf("%d",&num);
		for(i=0;i<num;i++)
			scanf("%d%d",&dat[i].a,&dat[i].b);
		qsort(dat,num,sizeof(node),cmpnode);
		idx=0;
		lis[0].next = 0;
		lis[0].b = dat[0].b;
		for(i=1;i<num;i++)
		{
			pcur=lis;b=dat[i].b;
			while(pcur && b<pcur->b)
			{ppre=pcur;pcur=pcur->next;}
			if(pcur)
				pcur->b=b;
			else
			{
				pcur = &lis[++idx];
				pcur->next=0;
				pcur->b=b;
				ppre->next=pcur;
			}
		}
		printf("%d\n",idx+1);
	}
	return 0;
}