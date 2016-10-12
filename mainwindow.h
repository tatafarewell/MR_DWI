#ifndef MAINWINDOW_H
#define MAINWINDOW_H

//// include ctk
//#include <ctkDICOMDatabase.h>
//#include <ctkDICOMIndexer.h>
//#include <ctkFileDialog.h>
////

//include QT
#include <QMainWindow>
//#include <QHash>
//#include <QString>
//#include <QStringList>
//#include <QProgressDialog>
//#include <QVariant>
//#include <QWidget>
//#include <QMessageBox>
//#include <qdebug.h>

namespace Ui {
class MainWindow;
}

class Ui_MainWindow;
class DicomModule;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

	/**
	* \brief CreateQtPartControl(QWidget *parent) sets the view objects from ui_QmitkDicomExternalDataWidgetControls.h.
	*
	* \param parent is a pointer to the parent widget
	*/
	//virtual void CreateQtPartControl(QWidget *parent);

	/**
	* \brief Initializes the widget. This method has to be called before widget can start.
	*/
	void Initialize();

signals:


	public slots :

	protected slots :
	void onStartdicom();
	void onDicomLoaded(bool);

private:


	Ui::MainWindow *ui;
	DicomModule * DicomUI;
};

#endif // MAINWINDOW_H
