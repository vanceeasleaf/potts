#include <iostream>
#include "potts.h"
#include "pottsSqure.h"
using namespace std;

int main()
{
    pottsSqure a;
    //Potts a;
    cout<<"������ Lx Ly q ��ʼ����(0=spinһ��,1=spin���:"<<endl;
    cin>>a.Lx>>a.Ly>>a.Nq>>a.init;
cout<<"�ӵ�nstart����ʼͳ�� ÿ��ִ��MCS�� ����(0=metropolis,1=sw,2=wolff):"<<endl;
cin>>a.nstart>>a.MCS>>a.method;
 cout<<"��ʼ�¶� �¶ȼ�� �����¶�:"<<endl;
cin>>a.t_start>>a.t_int>>a.t_end;
a.run();
    return 0;
}
