//
// Created by Арсений Юрченко on 25.10.2022.
//

#ifndef COSMIC_P4VECTOR_HH
#define COSMIC_P4VECTOR_HH


#include <iostream>
#include <math.h>
//using namespace std;
class p4vector {
private:
    double venergy, vx, vy, vz, vtheta, vphi, vmm;
public:
    p4vector(double senergy=.0, double sx=.0, double sy=.0, double sz=.0, double stheta=.0, double sphi=.0, double smm=.0) :
    venergy(senergy), vx(sx), vy(sy), vz(sz), vtheta(stheta), vphi(sphi), vmm(smm) {};
    double& e() { return venergy; }
    double& mm() { return vmm; }
    double& x() { return vx; }
    double& y() { return vy; }
    double& z() { return vz; }
    double& theta() { return vtheta; }
    double& phi() { return vphi; }
    double p() { return sqrt(vx*vx+vy*vy+vz*vz); }
    double mmr() { return venergy*venergy-vx*vx-vy*vy-vz*vz; }
    double thetar() { return acos(vz/(sqrt(vx*vx+vy*vy+vz*vz))); }
    double phir()
    {
        return atan2(vy, vx);
/*		if (vy/p()/sin(vtet)>.0)
		    {
		    if (vx/p()/sin(vtet)>.0) return asin(vy/p()/sin(vtet));
		    if (vx/p()/sin(vtet)<.0) return acos(vx/p()/sin(vtet));
		    }
		if (vy/p()/sin(vtet)<.0)
		    {
		    if (vx/p()/sin(vtet)>.0) return asin(vy/p()/sin(vtet));
		    if (vx/p()/sin(vtet)<.0) return -1.*acos(vx/p()/sin(vtet));
		    }*/
    }
    p4vector operator + (const p4vector a) { return p4vector(venergy+a.venergy, vx+a.vx, vy+a.vy, vz+a.vz); }
    p4vector operator - (const p4vector a) { return p4vector(venergy-a.venergy, vx-a.vx, vy-a.vy, vz-a.vz); }


};

//int main()
//{
//p4vec p1;
//p1.x()=p1.y()=p1.z()=1.;
//p1.tet()=p1.tetr();
//}


#endif //COSMIC_P4VECTOR_HH
