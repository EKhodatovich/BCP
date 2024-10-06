#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtCharts>
#include <QChartView>

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
};
#endif // MAINWINDOW_H


// class ChartView : public QChartView {
//     qreal mFactor=1.0;
// protected:
//     void wheelEvent(QWheelEvent *event) Q_DECL_OVERRIDE
//     {
//         chart()->zoomReset();

//         mFactor *= event->angleDelta().y() > 0 ? 0.5 : 2;

//         QRectF rect = chart()->plotArea();
//         QPointF c = chart()->plotArea().center();
//         rect.setWidth(mFactor*rect.width());
//         rect.moveCenter(c);
//         chart()->zoomIn(rect);

//         QChartView::wheelEvent(event);
//     }
// };
