#include "mainhub.h"
#include "ui_mainhub.h"


MainHub::MainHub(QWidget *parent) :
    QMainWindow(parent),
	ui(new Ui::MainWindow)
{
	ui->setupUi(this);
}

MainHub::~MainHub()
{
    delete ui;
}
