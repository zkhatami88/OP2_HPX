
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <fstream>
#include <cmath>
#include <array>
#include <hpx/hpx_init.hpp>
#include <hpx/runtime/actions/plain_action.hpp>
#include <hpx/hpx.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/parallel_algorithm.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/lcos/future_wait.hpp>
#include <hpx/lcos/wait_each.hpp>
#include <hpx/util/unwrapped.hpp>
#include <hpx/lcos/local/spinlock.hpp>
#include <hpx/exception.hpp>

#include <boost/timer.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>

#include <mutex>
#include <boost/thread.hpp>

#include "op_seq.h"

using hpx::naming::id_type;
using hpx::naming::invalid_id;

using hpx::lcos::future;
using hpx::lcos::wait;
using hpx::lcos::wait_each;
using hpx::async;

using hpx::util::high_resolution_timer;

using boost::program_options::variables_map;
using boost::program_options::options_description;
using boost::program_options::value;

using hpx::init;
using hpx::finalize;
using hpx::find_here;

using hpx::cout;
using hpx::flush;

using std::vector;
using namespace std;
typedef basic_ofstream<char> ofstream;
typedef unsigned int coord_t;



//Parameters for computing force and well seperated cells for each node
#define th 20

//Data for each node
struct Body {

public:
    int ID1,parent;
    double m1;
    std::vector<double> r1, v1, force; };
std::vector<Body> body;

void initialize_body(int N, Body body1){
    vector<double> A4(3,0);

    body1.ID1=0; body1.m1=0; body1.parent=0;
    body1.r1=A4; body1.v1=A4; body1.force=A4;

    for(int i=0; i<N; ++i)
        body.push_back(body1);}

struct Cell {

public:
    mutable hpx::lcos::local::spinlock mtx;
    int ID2, parent2,NumNodes=0,level,Ncell;
    vector<int> members, child,scell, list_cell1, list_cell2, neighbors;
    double m2;
    vector<double> r2, rd, boundary;


    void apply_changes(int);
    bool InCube(int, int );
    void det_boundary_subcube(int );

};
std::vector<Cell> cell(100000);

void* op_malloc (size_t size) {
  #ifdef __INTEL_COMPILER
    //return _mm_malloc(size, OP2_ALIGNMENT);
    return memalign(OP2_ALIGNMENT,size);
  #else
    return malloc(size);
  #endif
}

void initialize_cell(int N){
    vector<int> A2(6,-1); vector<double> A4(3,0),A5(6,0);

    cell[0].ID2=0; cell[0].parent2=0; cell[0].NumNodes=N; cell[0].level=0; cell[0].Ncell=0;
    cell[0].neighbors=A2; cell[0].r2=A4; cell[0].rd=A4; cell[0].boundary=A5;
    for (int i = 0; i < N; ++i)
        cell[0].members.push_back(i);}

//Read input from txt file
void read_Input(){
    int N; string temp,temp1;
    fstream textfile;
    string::size_type sz;
    textfile.open("ex10000.txt");
    textfile>>temp;
    N = atoi(temp.c_str());
    for (int i=0; i<N; ++i){
        textfile>>temp>>temp1;
        body[i].ID1=atoi(temp.c_str()); body[i].m1=stof(temp1,&sz);
        for (int j=2; j<5; ++j){
            textfile>>temp; body[i].r1[j-2]=stof(temp,&sz);}
        for (int j=5; j<8; ++j){
            textfile>>temp; body[i].v1[j-5]=stof(temp,&sz);}
        body[i].parent=0;}}

//Detremining max and min in the boundary cells
void Cell::apply_changes(int n){
    for(int j=0; j<cell[n].NumNodes; ++j){
        int i=cell[n].members[j];
        if (body[i].r1[0]<cell[n].boundary[0]) cell[n].boundary[0]=body[i].r1[0];
        if (body[i].r1[0]>cell[n].boundary[1]) cell[n].boundary[1]=body[i].r1[0];
        if (body[i].r1[1]<cell[n].boundary[2]) cell[n].boundary[2]=body[i].r1[1];
        if (body[i].r1[1]>cell[n].boundary[3]) cell[n].boundary[3]=body[i].r1[1];
        if (body[i].r1[2]<cell[n].boundary[4]) cell[n].boundary[4]=body[i].r1[2];
        if (body[i].r1[2]>cell[n].boundary[5]) cell[n].boundary[5]=body[i].r1[2];}}

//Finding out if the node is in that subcubic
bool Cell::InCube(int node, int p){
    for (int i=0; i<3; ++i)
        if(body[node].r1[i]<cell[p].boundary[2*i] || body[node].r1[i]>cell[p].boundary[2*i+1])
            return false;
    return true;}

//Determining the boundry of the each subcube
void Cell::det_boundary_subcube(int n){
    double a1,a2,a3,b1,b2,b3,c1,c2,c3;
    int p=n*8;

    a1=cell[n].boundary[0]; a2=(cell[n].boundary[0]+cell[n].boundary[1])/2; a3=cell[n].boundary[1];
    b1=cell[n].boundary[2]; b2=(cell[n].boundary[2]+cell[n].boundary[3])/2; b3=cell[n].boundary[3];
    c1=cell[n].boundary[4]; c2=(cell[n].boundary[4]+cell[n].boundary[5])/2; c3=cell[n].boundary[5];

    std::vector<double> A1, A2, A3, A4, A5, A6, A7, A8; std::vector<int> B2(6,-1),mem(cell[n].NumNodes); std::vector<double> B1(3,0);
    A1.push_back(a1); A1.push_back(a2); A1.push_back(b1); A1.push_back(b2); A1.push_back(c1); A1.push_back(c2);
    A2.push_back(a1); A2.push_back(a2); A2.push_back(b1); A2.push_back(b2); A2.push_back(c2); A2.push_back(c3);
    A3.push_back(a1); A3.push_back(a2); A3.push_back(b2); A3.push_back(b3); A3.push_back(c2); A3.push_back(c3);
    A4.push_back(a1); A4.push_back(a2); A4.push_back(b2); A4.push_back(b3); A4.push_back(c2); A4.push_back(c3);
    A5.push_back(a2); A5.push_back(a3); A5.push_back(b1); A5.push_back(b2); A5.push_back(c1); A5.push_back(c2);
    A6.push_back(a2); A6.push_back(a3); A6.push_back(b1); A6.push_back(b2); A6.push_back(c2); A6.push_back(c3);
    A7.push_back(a2); A7.push_back(a3); A7.push_back(b2); A7.push_back(b3); A7.push_back(c1); A7.push_back(c2);
    A8.push_back(a2); A8.push_back(a3); A8.push_back(b2); A8.push_back(b3); A8.push_back(c2); A8.push_back(c3);

    for(int j=1; j<9; ++j){ cell[p+j].r2=B1; cell[p+j].rd=B1; cell[p+j].neighbors=B2; cell[p+j].members=mem; }

    cell[p+1].boundary=A1; cell[p+2].boundary=A2; cell[p+3].boundary=A3; cell[p+4].boundary=A4;
    cell[p+5].boundary=A5; cell[p+6].boundary=A6; cell[p+7].boundary=A7; cell[p+8].boundary=A8;}

//creating octree and determining its parameters

void root(int n){
    vector<float> CM(3,0); cell[n].m2=0;
    for (int j=0; j<cell[n].NumNodes; ++j){
        cell[n].m2=cell[n].m2+body[cell[n].members[j]].m1;
        for (int i=0; i<3; ++i)
            CM[i]=CM[i]+body[cell[n].members[j]].r1[i] * body[cell[n].members[j]].m1;}
    for (int i=0; i<3;++i){
        cell[n].r2[i]=CM[i]/cell[n].m2;
        cell[n].rd[i]=cell[n].boundary[2*i+1]-cell[n].boundary[2*i];}}


/////////////////////////////////////////////////////////////////// Parellel Octree

void A(int n) {

    int p = n * 8;  Cell Foo;

    for (int i = 1; i < 9; ++i) {

        cell[p + i].ID2 = p + i;
        cell[p + i].level = cell[n].level + 1;
        cell[p + i].parent2 = n;

        for (int j = 0; j < cell[n].NumNodes; ++j) {
            if (Foo.InCube(cell[n].members[j], p + i)) {
                cell[p + i].members[cell[p + i].NumNodes] = cell[n].members[j];
                cell[p + i].NumNodes = cell[p + i].NumNodes + 1;}}

        if (cell[p + i].NumNodes > 1) root(p + i);

        if (cell[p + i].NumNodes >= 1 && cell[p + i].NumNodes <= th) {
            for (int j = 0; j < cell[p + i].NumNodes; ++j) {
                body[cell[p + i].members[j]].parent = n;
                cell[n].child.emplace_back(cell[p + i].members[j]);}}

        if (cell[p + i].NumNodes > th) {
            cell[n].scell.emplace_back(p + i);
            Foo.det_boundary_subcube(p + i);}}
}


hpx::future<void> strc3(int n){

    std::vector<hpx::future<void>> futures;
    futures.reserve(9);

    A(n);

    for(int i=1; i<9; ++i){
        int p=n*8;
        if(cell[p + i].NumNodes>th)
            futures.push_back(async(strc3,p+i));}

    return hpx::when_all(futures);
}

///////////////////////////////////////////////////////////////////


int hpx_main(int argc, char** argv){

    int N;
 /*   string temp;
    fstream textfile;
    textfile.open("ex10000.txt");
    textfile >> temp;
    N = (int) atoi(temp.c_str());

    Body body1;
    initialize_body(N, body1);     //Initializing body and cell:
    read_Input();
    Cell Foo;
    initialize_cell(N);
    Foo.apply_changes(0);
    root(0);
    Foo.det_boundary_subcube(0);

    ///////////////////////////////////////////////////////////

    boost::timer timer;
    hpx::future<void> r;
    //hpx::evaluate_active_counters(true, "rate_1");
    r=strc3(0);
    r.wait();
    double elapsed_time = timer.elapsed();
    std::cout << " CPU TIME: " << elapsed_time << endl;

    ////////////////////////////////////////////////////////////
   // hpx::id_type here = hpx::find_here();
    //std::cout<<here<<",";

    //hpx::evaluate_active_counters(true, "rate_2");

    /*std::cout<<endl;
    for(int i=0; i<100; ++i)
        std::cout<<cell[i].NumNodes<<","<<cell[i].ID2<<"-";*/

    std::cout<<"Yayyyy";

    hpx::finalize();
    return 0;
};

int main(int argc, char** argv){
    return hpx::init(argc, argv);
}
