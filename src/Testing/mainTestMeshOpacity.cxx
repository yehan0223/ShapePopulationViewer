//***************************************************************************//
//           Test setMeshOpacity in the class ShapePopulationBase            //
//***************************************************************************//

#include <iostream>
#include <string>
#include <QApplication>
#include <QFileInfo>

#include "testMeshOpacity.h"

int main(int, char *argv[])
{
    TestShapePopulationBase testShapePopulationBase;
    bool test = testShapePopulationBase.testMeshOpacity( (std::string)argv[1] );

    if(!test) return 0;
    else return -1;
}
