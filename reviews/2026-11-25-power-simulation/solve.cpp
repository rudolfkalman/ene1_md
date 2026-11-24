// Sparse conductance graph: grounded Laplacian, IC(0) preconditioned CG.
// Binary input: int32 N,M,K,ground; M*(int32 a,b,double G); K*int32 terminal.
// Outputs K arrays of N double node potentials, each for a 1 A injection.
#include <bits/stdc++.h>
using namespace std;
struct Edge { int a,b; double g; };
int main(int argc,char**argv){
 if(argc<3)return 1;
 ifstream f(argv[1],ios::binary);int n,m,k,ground;
 f.read((char*)&n,4);f.read((char*)&m,4);f.read((char*)&k,4);f.read((char*)&ground,4);
 vector<Edge> edges(m);f.read((char*)edges.data(),16LL*m);vector<int> terminals(k);f.read((char*)terminals.data(),4*k);
 if(!f)throw runtime_error("bad input");
 vector<map<int,double>> rows(n);vector<double>d(n,0);
 for(auto e:edges){d[e.a]+=e.g;d[e.b]+=e.g;if(e.a!=ground&&e.b!=ground){rows[e.a][e.b]-=e.g;rows[e.b][e.a]-=e.g;}}
 d[ground]=1;
 vector<int>ptr(n+1),col;vector<double>val;
 vector<vector<pair<int,double>>>lo(n);
 for(int i=0;i<n;i++){ptr[i]=col.size();for(auto [j,v]:rows[i]){col.push_back(j);val.push_back(v);if(j<i)lo[i].push_back({j,v});}map<int,double>().swap(rows[i]);}ptr[n]=col.size();
 vector<double>ld(n);
 for(int i=0;i<n;i++){
  double ss=0;
  for(int u=0;u<(int)lo[i].size();u++){
   int j=lo[i][u].first;double v=lo[i][u].second;int a=0,b=0;
   while(a<u&&b<(int)lo[j].size()){int x=lo[i][a].first,y=lo[j][b].first;if(x==y){v-=lo[i][a].second*lo[j][b].second;a++;b++;}else if(x<y)a++;else b++;}
   lo[i][u].second=v/ld[j];ss+=lo[i][u].second*lo[i][u].second;
  }
  ld[i]=sqrt(max(d[i]-ss,d[i]*1e-8));
 }
 auto mat=[&](const vector<double>&x,vector<double>&y){for(int i=0;i<n;i++){double s=d[i]*x[i];for(int a=ptr[i];a<ptr[i+1];a++)s+=val[a]*x[col[a]];y[i]=s;}};
 auto prec=[&](const vector<double>&r,vector<double>&z){for(int i=0;i<n;i++){double s=r[i];for(auto [j,v]:lo[i])s-=v*z[j];z[i]=s/ld[i];}for(int i=n-1;i>=0;i--){z[i]/=ld[i];for(auto [j,v]:lo[i])z[j]-=v*z[i];}};
 auto dot=[&](const vector<double>&a,const vector<double>&b){double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;};
 vector<double> supplied,mass;
 if(argc>=4){supplied.resize((long long)n*k);ifstream rf(argv[3],ios::binary);rf.read((char*)supplied.data(),8LL*n*k);if(!rf)throw runtime_error("bad RHS");}
 if(argc>=5){mass.resize(n);ifstream cf(argv[4],ios::binary);cf.read((char*)mass.data(),8LL*n);if(!cf)throw runtime_error("bad capacity");}
 ofstream out(argv[2],ios::binary);vector<double>previous(n,0);
 for(int c=0;c<k;c++){
  vector<double>x(n,0),r(n,0),z(n),p(n),ap(n),rhs(n,0);
  if(supplied.empty())rhs[terminals[c]]=1;else for(int i=0;i<n;i++)rhs[i]=supplied[(long long)c*n+i]+(mass.empty()?0:mass[i]*previous[i]);
  r=rhs;prec(r,z);p=z;double rz=dot(r,z);int it=0;
  for(;it<20000;it++){mat(p,ap);double alpha=rz/dot(p,ap);for(int i=0;i<n;i++){x[i]+=alpha*p[i];r[i]-=alpha*ap[i];}if(sqrt(dot(r,r))<1e-9)break;prec(r,z);double nrz=dot(r,z);double beta=nrz/rz;for(int i=0;i<n;i++)p[i]=z[i]+beta*p[i];rz=nrz;}
  mat(x,ap);for(int i=0;i<n;i++)ap[i]-=rhs[i];double energy=0;for(auto e:edges)energy+=e.g*pow(x[e.a]-x[e.b],2);
  cerr<<"case "<<c<<" n="<<n<<" edges="<<m<<" iter="<<it<<" R="<<setprecision(12)<<x[terminals[c]]<<" energy="<<energy<<" residual="<<sqrt(dot(ap,ap))<<endl;
  if(it>=20000||!isfinite(energy)||sqrt(dot(ap,ap))>1e-6)return 2;
  out.write((char*)x.data(),n*8LL);
  previous=x;
 }
}
