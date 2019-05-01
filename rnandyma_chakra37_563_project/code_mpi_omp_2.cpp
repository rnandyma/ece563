#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <queue>
#include <list>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <boost/mpi/status.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#include <boost/mpi.hpp>
#include<bits/stdc++.h>
#include <boost/optional.hpp>
namespace mpi = boost::mpi;
using namespace std;
omp_lock_t l0;
omp_lock_t l1;
omp_lock_t l2;
omp_lock_t l3;
omp_lock_t l4;
omp_lock_t l5;
omp_lock_t l6;


queue <std::string > fileQ;// filenames to be read
queue <std::string >opQ;//list of op files

const int nReaders=2;
const int nMappers=3;
const int nReducers=3;
const int nWriters=2;
int arr_reducers_done [nReducers];
int totalFileCount=0;
int linesRead=0;
int readersDone=0;
int linesMapped=0;
int loopcount = 0;
int mapperDone=0;
int reducerDone=0;
int proc_rank;
int numP;
//int mapsReduced = 0;
typedef struct

{

string content;

int count;

} record;

typedef struct
  {
   int node_id;
   int red_id;
  }identity;

typedef struct

{

string word;

int count;

} mapWord;

typedef queue<record> Qreader;  // Queue to go as input to the Mapper threads

Qreader rdQ;

void readFile(string fname) // reading characters from file and counting total

{

std::ifstream file(fname.c_str(), ios::app);
std::string str;



while (std::getline(file, str))
{

if(!str.empty()){

  transform(str.begin(), str.end(), str.begin(), ::tolower);
	str.push_back(' ');
	record r1;
        r1.content = str;
 r1.count  = strlen(str.c_str()); ;

 omp_set_lock(&l1);
 (rdQ).push(r1); linesRead++;
  
   omp_unset_lock(&l1);  
}
}

}

void reader(mpi::communicator world, int n)
{
string fname;
while(true){
 
world.send(numP-1,(100+n),"fName");
world.recv(numP-1,(100+n),fname);
if(fname.compare("NULL")==0)
	break;
else
 readFile(fname);
      

}
#pragma omp atomic
readersDone++;

}



map<string,int> wordMaps[nMappers];
vector <string>word_list;



void create_ip_files()
 {
  for (int i=1;i<=19;i++)

{

std::string name="abc";

name=name+std::to_string(i);

name+=".txt";


 fileQ.push(name);

 

}

 }

/*void create_op_files()
 {
  for (int i=0;i<nWriters;i++)

{
for (int j=0;j< numP-1;j++){
std::string name="output";

name=name+std::to_string(i)+std::to_string(j);

name+=".txt";


 opQ.push(name);

 

}}

 }*/

void insertMapperWords(string word, map<string,int> &tMap){

if (ispunct(word[word.size()-1])) 
        { 
            word.erase(word.size()-1, 1); }// remove punctuation if needed need to write an entire loop 

if(!tMap.empty()){
    if(tMap.count(word)>0){
	//omp_set_lock(&l2);
      tMap.find(word)->second++;
	//omp_unset_lock(&l2);
    }else{
	//omp_set_lock(&l2);
     tMap.insert({word,1});
	//omp_unset_lock(&l2);
    }

}else{
  
  tMap.insert({word,1});
	
}

 }

void mapper()
{
 int pid;
 string word;
 string s;
 map<string,int> myMap;
 
while(linesMapped<linesRead|| readersDone<nReaders)
{
   
    if(rdQ.empty()){
     usleep(500);
    
    } else{
       
            record r1;
           omp_set_lock(&l1);
       	 if(rdQ.empty())
	{
		omp_unset_lock(&l1);
 		continue;
	}
            r1=rdQ.front();
            (rdQ).pop();
            
            omp_unset_lock(&l1);
       
        s = r1.content.c_str();
             
             stringstream iss(s);
             
         
        // #pragma omp critical
        {
          while(iss >> word)
              {    
           
           
        
           insertMapperWords(word, myMap);
           
                  
              }}
        
	        #pragma omp atomic
		linesMapped++;
    }
		
}
omp_set_lock(&l2);
  wordMaps[mapperDone]=myMap;
  mapperDone++;
	omp_unset_lock(&l2);

}

identity hashing_f(string word, int nReducer, int numP)
{
 
  int j;

  int sum = 5381;
  
  identity id;
  int red_id;
  int hash, node_id;
  for(j=0;j<word.length(); j++)
   hash = ((sum << 5)+sum)+tolower(word[j]);
 
  hash = hash%(nReducer*(numP-1));
  id.node_id = hash%(numP-1);
  id.red_id = hash%nReducer;
 
  
  return id;
}


map<string, int>reducer_arr[nReducers];
void insertReducer(string word, int count, int red_id)
{
if(!reducer_arr[red_id].empty() && red_id<nReducers){
    if(reducer_arr[red_id].count(word)>0){
//omp_set_lock(&l2);
      reducer_arr[red_id].find(word)->second = (reducer_arr[red_id].find(word)->second)+count;
//omp_unset_lock(&l2);
    }else{
//omp_set_lock(&l2);
      reducer_arr[red_id].insert({word,count});
//omp_unset_lock(&l2);
    }

}else{
  //map<string,int> newreducer;
  //newreducer.insert({word,count});
omp_set_lock(&l3);
   reducer_arr[red_id].insert({word,count});
omp_unset_lock(&l3);
}

}

map<int, map<string, int> >red_bucket[20];


void insertNodeBucket(string word, int count, int red_id, int node_id, int numP)
{ // modify according to the node id and add the respective map to respective node so node_id will be checked and it will be sent to correct location
if(red_bucket[node_id].empty() && red_id<nReducers && node_id<numP-1){

omp_set_lock(&l6);
if(red_bucket[node_id].empty()){
map<string, int>* mymap= new map<string, int>();
     mymap->insert({word, count});
      red_bucket[node_id].insert({red_id,*mymap});
omp_unset_lock(&l6); return;}
else{
omp_unset_lock(&l6); 
}

}
if(!red_bucket[node_id].empty() && red_id<nReducers && node_id<numP-1){
    if(red_bucket[node_id].count(red_id)>0){
      if(red_bucket[node_id].find(red_id)->second.count(word)>0)
   {
   //omp_set_lock(&l4);
	red_bucket[node_id].find(red_id)->second.find(word)->second +=count;
   //omp_unset_lock(&l4);
   }
      else
  {
  //omp_set_lock(&l4);
	red_bucket[node_id].find(red_id)->second.insert({word,count});
  //omp_unset_lock(&l4);
  }
    }else{
      //omp_set_lock(&l5);
     map<string, int>* mymap= new map<string, int>();
     mymap->insert({word, count});
      red_bucket[node_id].insert({red_id,*mymap});
      //omp_unset_lock(&l5);

    }

}


}

void reducer(int n, mpi::communicator world)
{
  int pid;
  int currentMap;
  int mapsReduced = 0;
  string word;
  string s;
 // status stat = new status();
  identity recv_id;
  mpi::request reqs;
   int comm_complete = 0;
   map<string, int> new_map_recv;
	//int src=((proc_rank+1)%numP);
 // reqs =  world.irecv(src, n , new_map_recv);
	
  
while(mapperDone < nMappers || mapsReduced < mapperDone|| comm_complete < numP-2)
{
	boost::optional <mpi::status> stat= world.iprobe(mpi::any_source,n);
	if(stat && comm_complete < numP-2)
  	{
	mpi::status st= *stat;
	if(st.source()!=proc_rank){
	world.recv(st.source(), n , new_map_recv);
	
	if(!new_map_recv.empty()){
   	auto recv_it = new_map_recv.begin();
	while(recv_it != new_map_recv.end())
	{
	  insertReducer(recv_it->first,recv_it->second,n);
          recv_it++;
	}
	}
   comm_complete++;  	
	
   }
  }
  if(mapperDone == 0 || mapsReduced == mapperDone)
    { 
     //usleep(500);
    }
  else
    {  currentMap = mapsReduced;
       pid = omp_get_thread_num();
       record r2;
       
        
        for(auto map_it = wordMaps[currentMap].begin(); map_it != wordMaps[currentMap].end(); map_it++)
        {
         //omp_set_lock(&l3);
          r2.content = map_it->first.c_str();
          r2.count = map_it->second;
	  //omp_set_lock(&l3);
	// wordMaps[currentMap].erase(r2.content);
          //omp_unset_lock(&l3);
          recv_id = hashing_f(r2.content, nReducers, numP);
	   
	  if(recv_id.node_id == proc_rank)
          {
	  if(recv_id.red_id == n)
          insertReducer(r2.content, r2.count, recv_id.red_id);
	  else
          continue;
	  }
	  else{
	    if(recv_id.red_id == n)
            insertNodeBucket(r2.content, r2.count, recv_id.red_id, recv_id.node_id, numP);
            // now we have node_bucket,we need to send and recv and add to the respective red_id;
            // call function on to read from recv buffer and populate their respective reducer
	    else
	    continue;
		}
        }
        mapsReduced++;
     
      
     if(mapsReduced == nMappers)
       { map<std::string, int> new_map_send;

       for(int j = 0; j<numP-1;j++)
	{
	new_map_send.clear();
	if(j!=proc_rank){// only send for other nodes
       auto it = red_bucket[j].find(n);
       if(it != red_bucket[j].end()){
	 //new_map_send.insert(it->second.begin(),it->second.end());
	new_map_send=it->second;
      //world.send(j, n, it->second);
	world.send(j, n, new_map_send);
	
	}else{
	std::string str="xxxxx";
	new_map_send.insert({str,0});
	world.send(j, n, new_map_send);
        
	}
       }
       } 
      }
      
    }
}
omp_set_lock(&l4);
arr_reducers_done[reducerDone]=n;reducerDone++;
omp_unset_lock(&l4);


}
//string name2 = "output_small.txt";
void writeToFile(map<string, int> mapWrite, int n)
{

std::string name="output_old";

name=name+std::to_string(n)+std::to_string(proc_rank);

name+=".txt";
//opQ.push(name);
 ofstream fout;
 //string name2 = "output.txt";
 fout.open(name.c_str(), ios::app);
 //fout.open(name2.c_str(), ios::app);
if(!fout){
printf("error File");
return ;
}
for(auto map_it = mapWrite.begin(); map_it != mapWrite.end(); map_it++)
        {
                 if(map_it->first != "xxxxx")
		fout<< map_it->first<<"  " << map_it->second<<endl;
	}

fout.close();
}


void writer(int n){

int reducerWritten =0;

while(reducerDone<nReducers || reducerWritten < reducerDone){
 if(reducerWritten==reducerDone){
	usleep(500);
}else{

for (int i=reducerWritten;reducerWritten<reducerDone;i++){
if(arr_reducers_done[i]%nWriters==n){
writeToFile(reducer_arr[i],n);
}
reducerWritten++;
}


}

}
}



int main(int argc, char* argv[] )

{
double time1, time2;
int provided;
mpi::environment env(argc, argv,mpi::threading::level::multiple, true);
mpi::communicator world;
proc_rank = world.rank();
numP = world.size();
if(proc_rank==numP-1){
	omp_set_num_threads(1);
	create_ip_files();
	totalFileCount=fileQ.size();
	int null_sent=0; string fname;
	while(null_sent<((numP-1)*nReaders))
	{	
	boost::optional <mpi::status> stat= world.iprobe(mpi::any_source,mpi::any_tag);
	if(stat){
		mpi::status st=*stat; string rec;
		if(st.source()!=numP-1 && st.tag()>=100){
		world.recv(st.source(),st.tag() , rec);
		if(fileQ.empty())
		{
			fname ="NULL";
		world.send(st.source(),st.tag() , fname);
		null_sent++;
		} else
		{
  			fname =fileQ.front();
			world.send(st.source(),st.tag() ,fname);
  			 fileQ.pop();
		}
	}
	}		
	}
}else
{
omp_set_num_threads(5);
omp_init_lock(&l0);// fileQ
omp_init_lock(&l1);//rdQ
omp_init_lock(&l2);//mapperQ
omp_init_lock(&l3);
omp_init_lock(&l4);
omp_init_lock(&l5);
omp_init_lock(&l6);
time1 = omp_get_wtime();
#pragma omp parallel

{

#pragma omp master

{

for(int i = 0; i<nReaders;i++)
{
     #pragma omp task
reader(world,i);
}

for(int i = 0; i<nMappers;i++)
{
    #pragma omp task
mapper();
}

for(int i = 0; i<nReducers;i++)
{
   #pragma omp task
reducer(i, world);
}

for(int i = 0; i<nWriters;i++)
{
#pragma omp task
writer(i);
}

#pragma omp taskwait

}



}
time2 = omp_get_wtime();
printf("the time is %f\n", time2-time1);
omp_destroy_lock(&l0);
	omp_destroy_lock(&l1);
	omp_destroy_lock(&l2);
	omp_destroy_lock(&l3);
	omp_destroy_lock(&l4);
	omp_destroy_lock(&l5);
	omp_destroy_lock(&l6);
	
}

world.barrier();
env.~environment();
return 0;
}

