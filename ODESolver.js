class ODESolver {
  constructor(nVars, derivs) {
    this.nVars = nVars;
    this.Y = Array(nVars);
    this.dYdt = Array(nVars);
    this.K1 = Array(nVars);
    this.K2 = Array(nVars);
    this.K3 = Array(nVars);
    this.K4 = Array(nVars);
    this.work = Array(nVars);
    this.t = 0.0;
    let i=0;
    for(i=0;i<nVars;i++) {
        this.Y[i] = 0.0;
        this.dYdt[i] = 0.0;
        this.K1[i] = 0.0;
        this.K2[i] = 0.0;
        this.K3[i] = 0.0;
        this.K4[i] = 0.0;
        this.work[i] = 0.0;
    }
    this.derivs = derivs;
  }
  
  setInitialConditions(Y) {
      this.Y = Y.slice();
  }
  
  
  eulerStep(h) {
     let t = this.t;
     this.derivs(t,this.Y,this.dYdt);
     for(let i=0;i<this.nVars;i++) {
         this.Y[i] += h*this.dYdt[i];
     }
     this.t += h;
  }
  
  rk4Step(h) {
     let t = this.t;
     this.derivs(t,this.Y,this.dYdt);
     for(let i=0;i<this.nVars;i++) {
         this.K1[i] = h*this.dYdt[i];
         this.work[i] = this.Y[i]+this.K1[i]/2;
     }
     this.derivs(t+h/2,this.work,this.dYdt);
     for(let i=0;i<this.nVars;i++) {
         this.K2[i] = h*this.dYdt[i];
         this.work[i] = this.Y[i]+this.K2[i]/2;
     }
     this.derivs(t+h/2,this.work,this.dYdt);
     for(let i=0;i<this.nVars;i++) {
         this.K3[i] = h*this.dYdt[i];
         this.work[i] = this.Y[i]+this.K3[i];
     }
     this.derivs(t+h,this.work,this.dYdt);
     for(let i=0;i<this.nVars;i++) {
         this.K4[i] = h*this.dYdt[i];
         this.Y[i] = this.Y[i]+this.K1[i]/6+this.K2[i]/3+this.K3[i]/3+this.K4[i]/6;
     }
     this.t += h;

  }
  
  
  
}