#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define GRAPH_SIZE 490000
#define BOX_SIZE 200
#define ATOM_SIZE 3000
#define CONVEXH_SIZE 300
#define DEGREE_SIZE 565

long box[BOX_SIZE][BOX_SIZE][BOX_SIZE];
int satisfy[GRAPH_SIZE][3];
int compids[GRAPH_SIZE];
float atoms[ATOM_SIZE][4];
float catoms[ATOM_SIZE][4];
int cvh[CONVEXH_SIZE][3];
float cvhcenter[CONVEXH_SIZE][3];
float mins[CONVEXH_SIZE][3];
float maxs[CONVEXH_SIZE][3];
long Graph[GRAPH_SIZE][DEGREE_SIZE];
long Queue[GRAPH_SIZE];
long head;
long tail;
long ironNeighbors[GRAPH_SIZE];
// ******************************************************************
// min-max function
// ******************************************************************
float min3(float a,float b,float c)
{
      float m;
      
      m = a;
      if(b<m)
             m = b;
      if(c<m)
             m = c;
      return m;
}
float max3(float a,float b,float c)
{
      float m;
      
      m = a;
      if(b>m)
             m = b;
      if(c>m)
             m = c;
      return m;
}
int min(int a, int b)
{
    if(a<=b)
            return a;
    return b;
}
int max(int a, int b)
{
    if(a>=b)
            return a;
    return b;
}
// ******************************************************************
// queue function
// ******************************************************************
int initQueue()
{
    head = 0;
    tail = 0;
    return 0;
}
int isQempty()
{
    if(head==tail)
                  return 1;
    return 0;
}
int enQueue(long x)
{
    Queue[tail] = x;
    ++tail;
    return tail;
}
long deQueue()
{
     long y;
     y = Queue[head];
     ++head;
     return y;
}
long QLength()
{
     return tail-head;
}
// ******************************************************************
// main function
// ******************************************************************
int main( int argc, char* argv[] )
{
    // input prefix
    int joker = 0;
    char file_prefix[100];
    char file_atom[200];
    char file_corrected[200];
    char file_convex[200];
    char file_help[200];
    char file_box[200];
    char file_graph[200];
    FILE *fp;
    int acnt;
    int cacnt;
    float x;
    float y;
    float z;
    float r;
    int id1;
    int id2;
    int id3;
    int cvcnt;
    float minx;
    float maxx;
    float miny;
    float maxy;
    float minz;
    float maxz;
    float minr;
    float maxr;
    float R;
    float BR;
    float ironx;
    float irony;
    float ironz;
    float temp1;
    float temp2;
    int sx;
    int sy;
    int sz;
    float probe = 1.4;
    float d;
    float mind;
    int mindid;
    int intersect;
    float centx;
    float centy;
    float centz;
    long i;
    long j;
    long k;
    long lbid;
    float dist;
    int mindeg;
    int maxdeg;
    float avgdeg;
    int delta;
    float BR2;
    int ic;
    int jc;
    int kc;
    int lx;
    long components;
    long marked;
    float v1;
    float v2;
    float v3;
    float u1;
    float u2;
    float u3;
    float n1;
    float n2;
    float n3;
    float x1;
    float y1;
    float z1;
    float xx;
    float yy;
    float zz;
    int intersectCnt;
    float t;
    float ax;
    float by;
    float cz;
    int bibi;
    int bibj;
    int bibk;
    long indx;
     
    // ******************************************************************
    // input argument   
    // ******************************************************************
    memset(file_prefix, '\0', sizeof(file_prefix));
    memset(file_atom, '\0', sizeof(file_atom));
    memset(file_corrected, '\0', sizeof(file_corrected));
    memset(file_convex, '\0', sizeof(file_convex));
    memset(file_help, '\0', sizeof(file_help));
    memset(file_box, '\0', sizeof(file_box));
    memset(file_graph, '\0', sizeof(file_graph));
    if(argc>1)
    {
        strcpy(file_atom,argv[1]);
        strcpy(file_corrected,argv[1]);
        strcpy(file_convex,argv[1]);
        strcpy(file_help,argv[1]);
        strcpy(file_box,argv[1]);
        strcpy(file_graph,argv[1]);
    }
    else
    {
        strcpy(file_atom,"");
        strcpy(file_corrected,"");
        strcpy(file_convex,"");
        strcpy(file_help,"");
        strcpy(file_box,"");
        strcpy(file_graph,"");
    }
    strcat(file_atom,"all_atoms.txt");
    strcat(file_corrected,"all_correctedAtoms.txt");
    strcat(file_convex,"convexhull.txt");
    strcat(file_help,"helpings.txt");
    strcat(file_box,"boxes.txt");
    strcat(file_graph,"graph.txt");
    // ******************************************************************
    // read files  
    // ******************************************************************
    printf("Start reading %s:\n",file_atom);
    fp = fopen(file_atom, "r");
    if(fp==NULL)
    {
        printf("error in reading %s!\n",file_atom);
        return 0;
    }
    acnt = 0;
    while(!feof(fp))
    {
        fscanf(fp,"%f %f %f %f",&x,&y,&z,&r);
        atoms[acnt][0] = x;
        atoms[acnt][1] = y;
        atoms[acnt][2] = z;
        atoms[acnt][3] = r;
        ++acnt;
    }
    fclose(fp);
    printf("End of file-atoms: total number of atoms is %d\n\n",acnt);

    printf("Start reading %s:\n",file_corrected);
    fp = fopen(file_corrected, "r");
    if(fp==NULL)
    {
        printf("error in reading %s!\n",file_corrected);
        return 0;
    }
    cacnt = 0;
    while(!feof(fp))
    {
        fscanf(fp,"%f %f %f %f",&x,&y,&z,&r);
        catoms[cacnt][0] = x;
        catoms[cacnt][1] = y;
        catoms[cacnt][2] = z;
        catoms[cacnt][3] = r;
        ++cacnt;
    }
    fclose(fp);
    printf("End of file-adjusted coordinates: total number of atoms is %d\n\n",cacnt);

    printf("Start reading %s:\n",file_convex);
    fp = fopen(file_convex, "r");
    if(fp==NULL)
    {
        printf("error in reading %s!\n",file_convex);
        return 0;
    }
    cvcnt = 0;
    while(!feof(fp))
    {
        fscanf(fp,"%d %d %d",&id1,&id2,&id3);
        cvh[cvcnt][0] = id1;
        cvh[cvcnt][1] = id2;
        cvh[cvcnt][2] = id3;
        ++cvcnt;
    }
    fclose(fp);
    printf("End of file-convex hull: total number of atoms is %d\n\n",cvcnt);
    
    printf("Reading helping parameters from %s:\n",file_help);
    fp = fopen(file_help, "r");
    if(fp==NULL)
    {
        printf("error in reading %s!\n",file_help);
        return 0;
    }
    fscanf(fp,"%f %f",&minx,&maxx);
    fscanf(fp,"%f %f",&miny,&maxy);
    fscanf(fp,"%f %f",&minz,&maxz);
    fscanf(fp,"%f %f",&minr,&maxr);
    fscanf(fp,"%f %f",&R,&temp1);
    fscanf(fp,"%f %f",&temp1,&temp2);
    ironx = (temp1*1.0+temp2)/2.0;
    fscanf(fp,"%f %f",&temp1,&temp2);
    irony = (temp1*1.0+temp2)/2.0;
    fscanf(fp,"%f %f",&temp1,&temp2);
    ironz = (temp1*1.0+temp2)/2.0;
    fscanf(fp,"%f %f",&BR,&temp1);
    printf("X:[%2.3f %2.3f] Y:[%2.3f %2.3f] Z:[%2.3f %2.3f] R:[%minr %maxr]\n",minx,maxx,miny,maxy,minz,maxz,minr,maxr);
    printf("Iron: %2.3f %2.3f %2.3f\n",ironx,irony,ironz);
    printf("Box size: %2.3f Adjusted ball radius: %2.3f\n",R,BR);
    temp1 = probe;
    if(temp1<maxr)
        temp1 = maxr;
    sx = (int)round((maxx-minx+2*temp1)/R);
    sy = (int)round((maxy-miny+2*temp1)/R+1);
    sz = (int)round((maxz-minz+2*temp1)/R+1);
    printf("Large box size: %d %d %d\n\n",sx,sy,sz);
    // ******************************************************************
    // determine graph node candidates
    // ******************************************************************
    printf("Pre-calculations ...\n");
    // center of protein
    centx = 0;
    centy = 0;
    centz = 0;
    for(i=0;i<cacnt;++i)
    {
        centx += catoms[i][0];
        centy += catoms[i][1];
        centz += catoms[i][2];
    }
    centx = centx/cacnt;
    centy = centy/cacnt;
    centz = centz/cacnt;
    // center of cvh
    for(i=0;i<cvcnt;++i)
    {
        id1 = cvh[i][0];
        id2 = cvh[i][1];
        id3 = cvh[i][2];
        //v = (x2-x1 y2-y1 z2-z1) = v1 v2 v3
        //u = (x3-x1 y3-y1 z3-z1) = u1 u2 u3
        v1 = catoms[id2][0] - catoms[id1][0];
        v2 = catoms[id2][1] - catoms[id1][1];
        v3 = catoms[id2][2] - catoms[id1][2];
        u1 = catoms[id3][0] - catoms[id1][0];
        u2 = catoms[id3][1] - catoms[id1][1];
        u3 = catoms[id3][2] - catoms[id1][2];
        //(v2u3-u2v3  ,  u1v3-u3v1  ,  v1u2-v2u1)= (n1 n2 n3)
        cvhcenter[i][0] = v2*u3-u2*v3;
        cvhcenter[i][1] = v3*u1-u3*v1;
        cvhcenter[i][2] = v1*u2-u1*v2;
        mins[i][0] = min3(catoms[id1][0],catoms[id2][0],catoms[id3][0]);
        mins[i][1] = min3(catoms[id1][1],catoms[id2][1],catoms[id3][1]);
        mins[i][2] = min3(catoms[id1][2],catoms[id2][2],catoms[id3][2]);
        maxs[i][0] = max3(catoms[id1][0],catoms[id2][0],catoms[id3][0]);
        maxs[i][1] = max3(catoms[id1][1],catoms[id2][1],catoms[id3][1]);
        maxs[i][2] = max3(catoms[id1][2],catoms[id2][2],catoms[id3][2]);
    }
    // check empty spaces
    printf("Find empty spaces ... %d x %d x %d ...\n",sx,sy,sz);
    lbid = 0;
    for(i=0;i<sx;++i)
    {
        printf("%2.2f-",(i*1.0)/sx);
        for(j=0;j<sy;++j)
        {
            for(k=0;k<sz;++k)
            {
                box[i][j][k] = -1;
                // check if is inside or outside
                intersectCnt = 0;
                for(id1=0;id1<cvcnt;++id1)
                {
                    //t = (n1x1+n2x2+n3x3-n1x0-n2y0-n3z0)/n1
                    n1 = cvhcenter[id1][0];
                    n2 = cvhcenter[id1][1];
                    n3 = cvhcenter[id1][2];
                    id2 = cvh[i][0];
                    x1 = catoms[id2][0];
                    y1 = catoms[id2][1];
                    z1 = catoms[id2][2];
                    ax = 1;
                    by = 0;
                    cz = 0;
                    //t = (n1*x1+n2*y1+n3*z1-n1*i-n2*j-n3*k)/n1;
                    t = (n1*x1+n2*y1+n3*z1-n1*i-n2*j-n3*k)/(ax*n1+by*n2+cz*n3);
                    xx = i+t*ax;
                    yy = j+t*by;
                    zz = k+t*cz;
                    if(xx>=mins[id1][0] && xx<=maxs[id1][0])
                    {
                        if(yy>=mins[id1][1] && yy<=maxs[id1][1])
                        {
                            if(zz>=mins[id1][2] && zz<=maxs[id1][2])
                            {
                                ++intersectCnt;
                            }
                        }
                    }
                }
                if(intersectCnt%2==1)
                {
                    // if inside, any intersection with any atom?
                    intersect = 0;
                    for(id1=0;id1<cacnt;++id1)
                    {
                        dist = (catoms[id1][0]-i)*(catoms[id1][0]-i);
                        dist += (catoms[id1][1]-j)*(catoms[id1][1]-j);
                        dist += (catoms[id1][2]-k)*(catoms[id1][2]-k);
                        dist = sqrt(dist);
                        if(dist<(catoms[id1][3]+BR))
                        {
                            intersect = 1;
                            break;
                        }
                    }
                    if(intersect==0)
                    {
                        satisfy[lbid][0] = i;
                        satisfy[lbid][1] = j;
                        satisfy[lbid][2] = k;
                        compids[lbid] = 0;
                        Graph[lbid][0] = 0;
                        box[i][j][k] = lbid;
                        ++lbid;
                    }
                }
            }
        }
    }
    printf("\n");
    printf("Total number of (possible overlapping) empty spheres is %d!\n\n",lbid);
    // ******************************************************************
    // build graph
    // ******************************************************************
    printf("Build graph ... %d x %d x %d ...\n",sx,sy,sz);
    delta = (int)(round(BR)+1);
    BR2 = BR*BR;
    for(i=0;i<sx;++i)
    {
        printf("%2.2f-",(i*1.0)/sx);
        for(j=0;j<sy;++j)
        {
            for(k=0;k<sz;++k)
            {
                if(box[i][j][k]>-1)
                {
                    id1 = box[i][j][k];
                    for(ic=max(0,i-delta);ic<min(i+delta,sx);++ic)
                    {
                        for(jc=max(0,j-delta);jc<min(j+delta,sy);++jc)
                        {
                            for(kc=max(0,k-delta);kc<min(k+delta,sz);++kc)
                            {
                                if(box[ic][jc][kc]>-1)
                                {
                                    id2 = box[ic][jc][kc];
                                    d = (i-ic)*(i-ic)+(j-jc)*(j-jc)+(k-kc)*(k-kc);
                                    if(d<=BR2)
                                    {
                                        lx = Graph[id1][0];
                                        ++lx;
                                        Graph[id1][lx] = id2;
                                        Graph[id1][0] = lx;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("\n");
    mindeg = 20000000;
    maxdeg = -1;
    avgdeg = 0.0;
    for(id1=0;id1<lbid;++id1)
    {
        if(Graph[id1][0]<mindeg)
            mindeg = Graph[id1][0];
        if(Graph[id1][0]>maxdeg)
            maxdeg = Graph[id1][0];
        avgdeg += Graph[id1][0];
    }
    avgdeg = avgdeg/lbid;
    printf("Min-degree: %d  Max-degree: %d  Avg-degree: %2.3f \n\n",mindeg,maxdeg,avgdeg);
    // ******************************************************************
    // components in graph
    // ******************************************************************
    printf("Finding components in graph ...%d vertices ...\n",lbid);
    initQueue();
    marked = 0;
    enQueue(0);
    components = 1;
    compids[0] = components;
    ++marked;
    while(marked<lbid)
    {
        if(marked%1000==0)
            printf("%2.2f-",(marked*1.0)/lbid);
        if(isQempty()==0)
        {
            i = deQueue();
            for(j=1;j<=Graph[i][0];++j)
            {
                k = Graph[i][j];
                if(compids[k]==0)
                {
                    enQueue(k);
                    compids[k] = components;
                    ++marked;
                }
            }
        }
        else
        {
            ++components;
            for(i=0;i<lbid;++i)
            {
                if(compids[i]==0)
                {
                    enQueue(i);
                    compids[i] = components;
                    ++marked;
                    break;
                }
            }
        }
    }
    printf("\n");
    printf("Total number of components is: %d\n",components);
    // ******************************************************************
    // check iron
    // ******************************************************************
    printf("Build graph ... %d x %d x %d ...\n",sx,sy,sz);
    ironNeighbors[0] = 0;
    for(id1=0;id1<lbid;++id1)
    {
        ic = satisfy[id1][0];
        jc = satisfy[id1][1];
        kc = satisfy[id1][2];
        d = (ironx-ic)*(ironx-ic)+(irony-jc)*(irony-jc)+(ironz-kc)*(ironz-kc);
        if(d<=BR2)
        {
            lx = ironNeighbors[0];
            ++lx;
            ironNeighbors[lx] = id1;
            ironNeighbors[0] = lx;
        }
    }
    printf("Iron has %d neighors! \n\n",ironNeighbors[0]);
    //return(0);
    // ******************************************************************
    // writing results
    // ******************************************************************
    printf("Start writing %s:\n",file_box);
    fp = fopen(file_box, "w");
    if(fp==NULL)
    {
        printf("error in reading %s!\n",file_box);
        return 0;
    }
    acnt = 0;
    for(i=0;i<lbid;++i)
    {
        fprintf(fp,"%d %d %d %d %d\n",i,satisfy[i][0],satisfy[i][1],satisfy[i][2],compids[i]);
    }
    fclose(fp);
    printf("End of file-boxes\n\n");
    
    /* printf("Start writing %s:\n",file_graph);
    fp = fopen(file_graph, "w");
    if(fp==NULL)
    {
        printf("error in reading %s!\n",file_graph);
        return 0;
    }
    acnt = 0;
    for(i=0;i<lbid;++i)
    {
        fprintf(fp,"%d %d",i,Graph[i][0]);
        for(j=1;j<=Graph[i][0];++j)
            fprintf(fp,"%d ",Graph[i][j]);
        fprintf(fp,"\n");
    }
    fclose(fp);
    printf("End of file-graph\n\n");
    */
    return 0;
}
