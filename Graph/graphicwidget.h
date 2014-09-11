#ifndef GRAPHICWIDGET_H
#define GRAPHICWIDGET_H

#include <QWidget>
#include <QPen>
#include <QColor>
#include <QPainter>
#include <vector>
#include "graphdata.h"

class GraphicWidget : public QWidget
{
    Q_OBJECT
private:
    int iWidthLine;
    QPen pen;
    QColor col;
    std::vector<GraphData> GraphicsData;
     int w, h;
private: //methods
    void DrawLabels(QPainter &painter);

public: //methods
    GraphicWidget(QWidget *parent = 0);
    void AddGraphic(const GraphData &gd);
    void setColor(const QColor &color);
    void setLineWidth(int width);
    void paintEvent(QPaintEvent *e);
    void resizeEvent(QResizeEvent *e);
    ~GraphicWidget();
signals:
    
public slots:
    
};

#endif // GRAPHICWIDGET_H
