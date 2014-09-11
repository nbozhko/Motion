#include "graphicwidget.h"
#include "graphdata.h"

void GraphicWidget::DrawLabels(QPainter &painter)
{
    QFont sansFont("Helvetica [Cronyx]", 9);
    painter.setFont(sansFont);
    painter.setPen(QPen(Qt::black, 1));

    double aMin, aMax, fMin, fMax, coeffX, coeffY;


    for (std::vector<GraphData>::iterator it = GraphicsData.begin(); it != GraphicsData.end(); it++)
    {
        it->setWidthHeight(w, h);
    }
    coeffX = GraphicsData[0].getCoeffX();
    coeffY = GraphicsData[0].getCoeffY();

    aMin = GraphicsData[0].getArgMin();
    aMax = GraphicsData[0].getArgMax();
    fMin = GraphicsData[0].getValMin();
    fMax = GraphicsData[0].getValMax();

    double temp = aMin * coeffX;
    while(temp <= aMax * coeffX)
    {
        painter.drawLine(int( temp ), h, int( temp ), int(h - h*0.01));
        painter.drawText(int( temp ), int(h - h * 0.02), QString().sprintf("%.2f", temp / coeffX));
        temp += static_cast<double>(static_cast<int>((abs(aMax - aMin))  * coeffX / 10.0));
    }

    temp = fMin * coeffY;
    while(temp <= fMax * coeffY)
    {
       painter.drawLine(0, h - int( temp - fMin * coeffY), w * 0.01, h - int( temp - fMin * coeffY ));
       painter.drawText(w * 0.02, h - int( temp - fMin * coeffY), QString().sprintf("%.2f", temp / coeffY));
       temp += static_cast<double>(static_cast<int>((abs(aMax - aMin))  * coeffX / 10.0));
    }

}

void GraphicWidget::setLineWidth(int width)
{
    iWidthLine = width;
    pen.setWidth(iWidthLine);
}

void GraphicWidget::setColor(const QColor &color)
{
    col = color;
    pen.setColor(col);
}

GraphicWidget::GraphicWidget(QWidget *parent) :
    QWidget(parent)
{
    this->setStyleSheet("background-color: white");
    iWidthLine = 1;
    w = this->width();
    h = this->height();
}

void GraphicWidget::paintEvent(QPaintEvent *e)
{
    QPainter painterWidget;
    w = this->width();
    h = this->height();
    for (std::vector<GraphData>::iterator it = GraphicsData.begin(); it != GraphicsData.end(); it++)
    {
        it->setWidthHeight(w, h);
    }

    painterWidget.begin(this);
        painterWidget.setPen(pen);
        painterWidget.setRenderHint(QPainter::HighQualityAntialiasing, true);

        for (std::vector<GraphData>::iterator it = GraphicsData.begin(); it != GraphicsData.end(); it++)
        {
            for (unsigned i = 0; i < it->numPoints() - 1; i++)
            {
                painterWidget.drawLine(int(it->getArg(i)), int(it->getVal(i)),
                                       int(it->getArg(i + 1)), int(it->getVal(i + 1)));
            }
        }
        DrawLabels(painterWidget);
        painterWidget.drawRect(0, 0, int(w - 0.001 * w), int(h - h * 0.001));
    painterWidget.end();
}

void GraphicWidget::resizeEvent(QResizeEvent *e)
{
    this->setStyleSheet("background-color: white");
    w = this->width();
    h = this->height();
    for (std::vector<GraphData>::iterator it = GraphicsData.begin(); it != GraphicsData.end(); it++)
    {
        it->setWidthHeight(w, h);
    }
}

GraphicWidget::~GraphicWidget()
{
    iWidthLine = 0;
}

void GraphicWidget::AddGraphic(const GraphData &gd)
{
    GraphicsData.push_back(gd);
}
