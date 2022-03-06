// by Concyclics
//
//
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include <thread>
#include <chrono>
#include <time.h> 
#include <random>
#include <mutex>

//#define time_show
//#define performance_show

//#define mem_clear

#define envalueFunction

//#define back_edge_cnt

using namespace std;


const int thread_cnt=32;
const double eps=1e-5;//presicion error
const int loop_min=3;
const int loop_max=6;
const double low_R=0.9;
const double high_R=1.1;

int cut_cnt;

mt19937 rnd(time(NULL));

int64_t NowTimeMillis()
{
    int64_t timems = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    return timems;
}

struct node//def of a edge
{
    //int begin;
    unsigned int end;
    int value;
    unsigned long long time;
    
    //compare by time to use sort() and upper_bound
    bool operator <(const node &x) const;
    bool operator >(const node &x) const;
    bool operator <=(const node &x) const;
    bool operator >=(const node &x) const;
    bool operator ==(const node &x) const;
    
#ifdef performance_show
    void display();
#endif
};

unordered_map<long long,unsigned int> IDcheck;//account ID, key
unordered_map<unsigned int,long long> IDcheckback;//key account ID

int N,M;//number of accounts

vector<vector<node> >edge;//graph
vector<node> graph;
vector<unsigned  int> head;
vector<unsigned int> dfn_order;
vector<char> vis;
vector<bool> edge_lock;
int tMc[thread_cnt];


#ifdef back_edge_cnt
vector<node> back_edge[account_MAX];//graph
volatile bool back_edge_lock[account_MAX];
#endif

vector<vector<node> > answers[thread_cnt];//store all answers


//read buffers
vector<char> * input;
long long pos=0;
long long pos_mult[thread_cnt+1];

double MAX=0;

mutex account_finished;

//read Account from file and make keys from account.csv
void readAccount(char* path);

//read transfer from file and sort by time from transfer.csv
void readEdge(char* path);

void read_transfer_spilit();

//DFS thread
bool DFS(const int begin, int now, int t, const int thread_NO, node ans[]);

inline bool inLoop(const int x, const int t, node ans[]);

//use DFS from a point and search answers
void one_thread_DFS(int x);

//multiprocess sort
void one_thread_Sort(int x);


void one_thread_read(int p);

//multicore search
void findAns_multicore();

void sort_multicore();

void read_multicore();

//print all answer to output.csv
void printAns(char* path);

//fast read from input vector
inline bool read_from_buff_int64(long long &b);

inline void read_from_buff_double(int &b);

//fast read for multicore
inline bool read_from_buff_int64_mult(long long &b, long long &l, long long &r);

inline void read_from_buff_double_mult(int &b, long long &l, long long &r);

void dfn(int x)
{
    vis[x]=(char)1;
    dfn_order.push_back(x);
    int i,v;
    for (i=head[x];i<head[x+1];++i)
    {
        v=graph[i].end;
        if (!vis[v])
            dfn(v);
    }
}

bool envalueFvalue(int V)
{
    return V>=23600&&V<=1000000;
}

bool envalueFtime(unsigned long long T)
{
    return T<1358431033037ull&&T>1262660165959ull;
}


int main(int argc, char *argv[])
{
    

   char * account_path="./scale1/account.csv";
   char * transfer_path="./scale1/transfer.csv";
   char * output_path="output.csv";

    time_t start = NowTimeMillis();
    
    if(argc==4)
    {
        account_path=argv[1];
        transfer_path=argv[2];
        output_path=argv[3];
    }
    
    std::thread read_account_thread=std::thread(readAccount, account_path);
    
    std::thread read_transfer_thread=std::thread(readEdge, transfer_path);
    
    read_account_thread.join();
    
    read_transfer_thread.join();
    
#ifdef  performance_show
    int64_t copyT=NowTimeMillis();
#endif
    
    head.resize(N+1);
    for(int i=1;i<N;++i)
        head[i+1]=head[i]+edge[i].size();
    M=head[N];
    graph.resize(M+1);
    
    for (int i=0;i<N;++i)
    {
        copy(edge[i].begin(),edge[i].end(),graph.begin()+head[i]);
#ifdef mem_clear
        edge[i].shrink_to_fit();
#endif
    }

#ifdef performance_show
    int64_t copyT2=NowTimeMillis();
    cout<<"copy end "<<copyT2-copyT<<"ms\n";
#endif    
    dfn_order.reserve(N+1);
    vis.resize(N+1);
    for (int i=0;i<N;++i)
        if (!vis[i])
            dfn(i);
    
    
#ifdef time_show
    
    time_t read_end = NowTimeMillis();
    
    cout<<'\n'<<"read file end\nUsing time "<<(int)(read_end-start)<<"ms\n\n";
    
#endif
    
    
    findAns_multicore();
    
#ifdef time_show
    
    time_t find_end = NowTimeMillis();
    
    cout<<"find answer end\nUsing time "<<(int)(find_end-read_end)<<"ms\n\n";
    
#endif
    
    printAns(output_path);
    
#ifdef time_show
    
    time_t output_end = NowTimeMillis();
    
    cout<<"output end\nUsing time "<<(int)(output_end-find_end)<<"ms\n\n";
    
    cout<<"total using time "<<(int)(output_end-start)<<"ms\n\n";
    
    
    cout<<cut_cnt;
#endif
    
}


void readAccount(char* path)
{
    
    account_finished.lock();
    
    ifstream fin(path, std::ios::binary);
    
    vector<char> buf(fin.seekg(0, std::ios::end).tellg());
    fin.seekg(0, std::ios::beg).read(&buf[0], static_cast<std::streamsize>(buf.size()));
    
    input=&buf;
    
    fin.close();
    
    if(*(input->end()-1)!='\n')input->emplace_back('\n');
    
    int cnt=0;
    long long tmp;
    
//    cout<<"buff finished \n";
    
    pos=0;
    
    while(read_from_buff_int64(tmp))
    {
        IDcheckback[cnt]=tmp;
        IDcheck[tmp]=cnt;//make key
        ++cnt;
        //printf("reading account : %d \r", cnt);
    }
    
    //printf("reading account : %d \n", cnt);
    
    //cout<<'\n';
    
    account_finished.unlock();
#ifdef mem_clear
    input->shrink_to_fit();
#endif
    N=cnt;
}


void readEdge(char* path)
{
#ifdef time_show
    
    long long fro=NowTimeMillis();
    
#endif
    
    ifstream fin(path, std::ios::binary);
    
    vector<char> buf(fin.seekg(0, std::ios::end).tellg());
    fin.seekg(0, std::ios::beg).read(&buf[0], static_cast<std::streamsize>(buf.size()));
    
    fin.close();
    
    input=&buf;
    
    if(*(input->end()-1)!='\n')input->emplace_back('\n');
    
#ifdef time_show
    
    long long from=NowTimeMillis();
    
    cout<<"process read transfer use "<<from-fro<<"ms\n";
    
#endif
    
    read_transfer_spilit();
    
    account_finished.lock();
    
    edge.resize(N);
    
    edge_lock.resize(N);
    
    read_multicore();
    
    account_finished.unlock();
#ifdef mem_clear
    input->shrink_to_fit();
#endif
    
#ifdef time_show
    
    long long kkk=NowTimeMillis();
    
    cout<<"process fast read use "<<kkk-from<<"ms\n";
    
#endif
    
    sort_multicore();
    
#ifdef time_show
    
    long long to=NowTimeMillis();
    
    cout<<"process sort use "<<to-kkk<<"ms\n";
    
#endif
    
    int cnt=0;
    
    //printf("reading transfer : %d \r",cnt);
    
}

void read_transfer_spilit()
{
    long long total_size=input->size();
    
    long long each_size=(total_size)/thread_cnt;
    
    pos_mult[0]=0;
    
    pos_mult[thread_cnt]=total_size;
    
    for(int i=1;i<thread_cnt;++i)
    {
        pos_mult[i]=pos_mult[i-1]+each_size;
    }
    
    for(int i=thread_cnt-1;i>0;--i)
    {
        while((*input)[pos_mult[i]]!='\n')++pos_mult[i];
    }
    
    return;
    
}


bool DFS(const int begin, int now, int t, const int thread_NO, node ans[])
{
    
    if(t>=loop_min&&ans[t-1].end==begin)
    {
        answers[thread_NO].emplace_back(ans, ans+t);
        return true;
    }
    
    
    if(t==loop_max)
    {
        return false;
    }

    bool ret=false;

    for(int i=head[now];i<head[now+1];++i)//get the first edge timestamp>now
    {
        if(graph[i].time<=ans[t-1].time)break;
        
        if(abs(graph[i].value-ans[t-1].value)*10<=ans[t-1].value&&!inLoop(graph[i].end,t,ans))//check value in [0.9,1.1] and point do not visited
        {
            ans[t]=graph[i];
            ret|=DFS(begin,graph[i].end,t+1,thread_NO,ans);
            //ans.pop_back();
        }
    }
    
    return ret;
}

bool inLoop(const int x, const int t, node ans[])
{
    for(int i=0;i<t;i++)
    {
        if(ans[i].end==x)return 1;
    }
    return 0;
}


void one_thread_DFS(int x)
{
    //vector<node> ans;
    for(int I=x;I<N;I+=thread_cnt)
    {
        unsigned int i=dfn_order[I];
        bool ret=false;
#ifdef back_edge_cnt
        unsigned long long MAXtime=0;;
        int MAXvalue=0, MINvalue=1<<30;
        
        for(auto &X:back_edge[i])
        {
            MAXvalue=max(MAXvalue, X.value);
            MINvalue=min(MINvalue, X.value);
            MAXtime=max(MAXtime, X.time);
        }
#endif
        for(int j=head[i];j<head[i+1];++j)
        {
#ifdef back_edge_cnt
            if(graph[j].value*1.3<MINvalue||graph[j].value*0.7>MAXvalue||graph[j].time+1000>MAXtime)
            {
                cut_cnt++;
                continue;
            }
#endif
            
#ifdef envalueFunction
            if(!envalueFvalue(graph[j].value)||!envalueFtime(graph[j].time))
            {
                cut_cnt++;
                continue;
            }
#endif
            //if(ret&&rnd()&1)continue;
            node ans[loop_max];
            ans[0]=graph[j];
            //ans.emplace_back(T);
            ret|=DFS(i,graph[j].end,1,x,ans);
            //ans.pop_back();
        }
        
        //if(!x)printf("working...  %.2lf%%\r", (i+1)*100.0/N);
    }
    return;
}

void one_thread_Sort(int x)
{
    
    for(int i=x;i<N;i+=thread_cnt)
    {
        sort(edge[i].begin(),edge[i].end(),greater<node>());
    }
    return;
}

void one_thread_read(int p)
{
    
    long long l=pos_mult[p];
    long long r=pos_mult[p+1];
    
    //printf("thread %d l %lld r %lld \n",p,l,r);
    
    long long x;
    int begin,end;
    int X;
    node Etmp;
    
    while(read_from_buff_int64_mult(x,l,r))
    {
        begin=IDcheck[x];
        //Etmp.begin=begin;
        
        read_from_buff_int64_mult(x,l,r);
        end=IDcheck[x];
        Etmp.end=end;
        
        read_from_buff_int64_mult(x,l,r);
        Etmp.time=x;
        
        read_from_buff_double_mult(X,l,r);
        Etmp.value=X;
        

        while(edge_lock[begin]);
        edge_lock[begin]=true;
        edge[begin].emplace_back(Etmp);
        edge_lock[begin]=false;
        
        
#ifdef back_edge_cnt
        
        while(back_edge_lock[end]);
        back_edge_lock[end]=true;
        back_edge[end].emplace_back(Etmp);
        back_edge_lock[end]=false;
        
#endif        
        //printf("reading transfer : %d \r",cnt);
        //Etmp.display();
    }
    //printf("thread %d end\n",p);
    return;
}


void findAns_multicore()
{
    
    thread ID[thread_cnt];
    
    for(int i=thread_cnt-1;i>=0;--i)
    {
        ID[i]=thread(one_thread_DFS,i);
    }
        
    for(int i=thread_cnt-1;i>=0;--i)
    {
        ID[i].join();
    }
    
    //printf("working...  %.2lf%%\n",100.0);
    
    //cout<<'\n';
}

void sort_multicore()
{
    thread ID[thread_cnt];
    
    for(int i=thread_cnt-1;i>=0;--i)
    {
        ID[i]=thread(one_thread_Sort,i);
    }
    
    for(int i=thread_cnt-1;i>=0;--i)
    {
        ID[i].join();
    }


}

void read_multicore()
{
    thread ID[thread_cnt];
    
    for(int i=thread_cnt-1;i>=0;--i)
    {
        ID[i]=thread(one_thread_read,i);
    }
    
    for(int i=thread_cnt-1;i>=0;--i)
    {
        ID[i].join();
    }

    
}


void printAns(char* path)
{
    
    std::ofstream output(path);
    
    int total=0;
    
    
    for(int cnt=0;cnt<thread_cnt;++cnt)
    {
        for(auto &x:answers[cnt])
        {
            output<<'('<<IDcheckback[(x.end()-1)->end]<<')';
            //Rr=x[0].value;
            //cout<<(x.end()-1)->end<<endl;
            for(auto &p:x)
            {
                output<<"-["<<p.time<<',';
                output<<std::fixed<<std::setprecision(2)<<p.value/1000.0;
                output<<"]->("<<IDcheckback[p.end]<<')';
            }
            output<<'\n';
            total++;
            //printf("output answer : %d \r",total);
        }
    }
    
#ifdef performance_show
    
    printf("output answer : %d \n",total);
    
    
#endif
    
    output.close();
    
}

inline bool read_from_buff_int64(long long &b)
{
    b=0;
    char a=(*input)[pos];
    while ( (a>'9' || a<'0') && pos<input->size())
        a=(*input)[++pos];
    while ( (a>='0' && a<='9') && pos<input->size())
    {
        b=b*10+a-'0';
        a=(*input)[++pos];
    }
    return pos<input->size();
}

inline bool read_from_buff_int64_mult(long long &b, long long &l, long long &r)
{
    b=0;
    char a=(*input)[l];
    while ( (a>'9' || a<'0') && l<r)
        a=(*input)[++l];
    while ( (a>='0' && a<='9') && l<r)
    {
        b=b*10+a-'0';
        a=(*input)[++l];
    }
    return l<r;
}

/*
inline void read_from_buff_double(double &b)
{
b=0;
char a=(*input)[pos];
int c=0;
while ( (a>'9' || a<'0') && pos<input->size())
{
a=(*input)[++pos];
}

while ((a>='0' && a<='9')&&pos<input->size())
{
b=b*10+(a-'0');
a=(*input)[++pos];
}
if(a=='.')
{
a=(*input)[++pos];
while ((a>='0' && a<='9')&&pos<input->size())
{
b=b*10+(a-'0');
a=(*input)[++pos];
++c;
}
}
while (c--)b*=0.1;
return;
}
*/

inline void read_from_buff_double(int &b)
{
    b=0;
    char a=(*input)[pos];
    int c=0;
    while ( (a>'9' || a<'0') && pos<input->size())
    {
        a=(*input)[++pos];
    }

    while ((a>='0' && a<='9')&&pos<input->size())
    {
        b=b*10+(a-'0');
        a=(*input)[++pos];
    }
    if(a=='.')
    {
        a=(*input)[++pos];
        while ((a>='0' && a<='9')&&pos<input->size())
        {
            b=b*10+(a-'0');
            a=(*input)[++pos];
            ++c;
        }
    }
    while (++c<=3)b*=10;
    return;
}

/*
inline void read_from_buff_double_mult(double &b, long long &l, long long &r)
{
b=0;
char a=(*input)[l];
int c=0;
while ( (a>'9' || a<'0') && l<r)
{
a=(*input)[++l];
}

while ((a>='0' && a<='9')&&l<r)
{
b=b*10+(a-'0');
a=(*input)[++l];
}
if(a=='.')
{
a=(*input)[++l];
while ((a>='0' && a<='9')&&l<r)
{
b=b*10+(a-'0');
a=(*input)[++l];
++c;
}
}
while (c--)b*=0.1;
return;
}
*/

inline void read_from_buff_double_mult(int &b, long long &l, long long &r)
{
    b=0;
    char a=(*input)[l];
    int c=0;
    while ( (a>'9' || a<'0') && l<r)
    {
        a=(*input)[++l];
    }
    
    while ((a>='0' && a<='9')&&l<r)
    {
        b=b*10+(a-'0');
        a=(*input)[++l];
    }
    if(a=='.')
    {
        a=(*input)[++l];
        while ((a>='0' && a<='9')&&l<r)
        {
            b=b*10+(a-'0');
            a=(*input)[++l];
            ++c;
        }
    }
    while (++c<=3)b*=10;
    return;
}


bool node::operator <(const node &x) const
{
    return time<x.time;
}

bool node::operator >(const node &x) const
{
    return time>x.time;
}

bool node::operator <=(const node &x) const
{
    return time<=x.time;
}

bool node::operator >=(const node &x) const
{
    return time>=x.time;
}

bool node::operator ==(const node &x) const
{
    return time==x.time;
}

#ifdef performance_show
void node::display()
{
    //cout<<"from:"<<begin<<" to:"<<end<<" intime:"<<time<<" byvalue:"<<value<<"\n";
}
#endif
