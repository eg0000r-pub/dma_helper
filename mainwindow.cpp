#include <cmath>
#include <iostream>

#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "fsolve.h"

#ifndef M_PI
#define M_PI (4.0*atan(1.0))
#endif

#define E_C (1.602e-19)

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->tSet->setValue(298.0);
    ui->pSet->setValue(1.0);
    ui->dpSet->setValue(0.15);
    ui->nSet->setValue(1);
    ui->qshSet->setValue(6.5);
    ui->lSet->setValue(0.44369);
    ui->r2Set->setValue(0.01961);
    ui->r1Set->setValue(0.00937);
    calcV();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::calcV()
{
    double t, tr, p, pr, dp, n, q, l, r1, r2, cc, zp, mfp, mu, v;
    t = ui->tSet->value();
    tr = t/298.0;
    p = ui->pSet->value();
    pr = p/1.0;
    dp = ui->dpSet->value()*1.0e-6;
    n = double(ui->nSet->value());
    q = ui->qshSet->value();
    l = ui->lSet->value();
    r1 = ui->r1Set->value();
    r2 = ui->r2Set->value();
    mfp = tr*0.0000000651/pr;
    cc = 1+(2*mfp/dp)*(1.257+0.4*exp(-1.1*dp/(2*mfp)));
    mu = t*0.0000000464+0.00000452;
    zp = n*E_C*cc/(mu*3*M_PI*dp);
    v = -((q/60000)*(log(r1/r2))/(2*M_PI*l*zp));
    ui->vSet->setValue(v);
}

void MainWindow::on_calcBtn_released()
{
    calcV();
}


void MainWindow::on_solveBtn_released()
{
    double t, tr, p, pr, dp0, n, q, l, r1, r2, mfp, mu;
    t = ui->tSet->value();
    tr = t/298.0;
    p = ui->pSet->value();
    pr = p/1.0;
    dp0 = ui->dpSet->value()*1.0e-6;
    n = double(ui->nSet->value());
    q = ui->qshSet->value();
    l = ui->lSet->value();
    r1 = ui->r1Set->value();
    r2 = ui->r2Set->value();
    mfp = tr*0.0000000651/pr;
    mu = t*0.0000000464+0.00000452;
    auto f = [this,n,mfp,mu,r1,r2,l,q](double dp) -> double {
        double cc, zp, v;
        cc = 1+(2*mfp/dp)*(1.257+0.4*exp(-1.1*dp/(2*mfp)));
        zp = n*E_C*cc/(mu*3*M_PI*dp);
        v = -((q/60000)*(log(r1/r2))/(2*M_PI*l*zp));
        return v-ui->vTarget->value();
    };
    double dp = fsolve(f, dp0);
    ui->dpSet->setValue(dp*1.0e6);
    calcV();
    //ui->dpSet->setValue(0.240);
    //std::cout << calcV()-ui->vTarget->value() << std::endl;
}

