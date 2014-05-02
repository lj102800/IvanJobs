#include <stdio.h>

int main()
{
	int a[5000],b[5000];
	int num,i,j,idx,tmp,tail,count;
	scanf("%d",&count);
	while(count>0)
	{
		count--;
		scanf("%d",&num);
		for(i=0;i<num;i++)
			scanf("%d%d",&a[i],&b[i]);
		for(i=0;i<num;i++)
		{
			idx=i;
			for(j=i+1;j<num;j++)
				if(a[idx]>a[j]||(a[idx]==a[j]&&b[idx]>b[j]))
					idx=j;
			tmp=a[i];
			a[i]=a[idx];
			a[idx]=tmp;
			tmp=b[i];
			b[i]=b[idx];
			b[idx]=tmp;
		}
		tail=0;
		for(i=1;i<num;i++)
		{
			j=tail;
			while(j>=0&&b[i]<b[j])
				j--;
			if(j!=-1)
				b[j]=b[i];
			else
			{
				tail++;
				tmp=b[i];
				for(j=tail;j>0;j--)
					b[j]=b[j-1];
				b[0]=tmp;
			}
		}
		printf("%d\n",tail+1);
	}
	return 0;
}