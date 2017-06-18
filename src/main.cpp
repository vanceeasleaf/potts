#include <iostream>
#include "potts.h"
#include "pottsSqure.h"
using namespace std;

int main()
{
	pottsSqure a;
	//Potts a;
	cout<<"PLEASE INPUT : Lx,Ly,q,init_condiction(0: spin parallel,1:spin random)"<<endl;
	cin>>a.Lx>>a.Ly>>a.Nq>>a.init;
	cout<<"PLEASE INPUT : nstart,MCS,method(0=metropolis,1=sw,2=wolff):"<<endl;
	cin>>a.nstart>>a.MCS>>a.method;
	 cout<<"PLEASE INPUT : t_start,t_interval,t_end"<<endl;
	cin>>a.t_start>>a.t_int>>a.t_end;
	a.run();
    return 0;
}
