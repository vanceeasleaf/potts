#include <iostream>
#include "potts.h"
#include "pottsSqure.h"
using namespace std;

int main()
{
    pottsSqure a;
    //Potts a;
    cout<<"请输入 Lx Ly q 初始条件(0=spin一致,1=spin随机:"<<endl;
    cin>>a.Lx>>a.Ly>>a.Nq>>a.init;
cout<<"从第nstart步开始统计 每次执行MCS步 方法(0=metropolis,1=sw,2=wolff):"<<endl;
cin>>a.nstart>>a.MCS>>a.method;
 cout<<"开始温度 温度间隔 结束温度:"<<endl;
cin>>a.t_start>>a.t_int>>a.t_end;
a.run();
    return 0;
}
