#!/usr/bin/env python
# Visualizing data from a molecular dynamics simulation
# Data provided by Daniel Spangberg
from vtk import *
import argparse, os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir))
#add parentdir into PYTHONPATH, where IO module can be found
import IO

version=vtkVersion().GetVTKMajorVersion()

def Read(filename):
    Dict=IO.LoadDict(filename)
    data=Dict["Points"]
    points = vtkPoints()
    color = vtkFloatArray()
    captions = []
    SubLat=set()
    for vec, sub in data:
        x, y = vec[0], vec[1]
        if len(vec)==3:
            z=vec[2]
        else:
            z=0
        points.InsertNextPoint(x, y, z)
        color.InsertNextValue(sub)
        captions.append([(x,y,z),str(tuple(coord))])
        SubLat.add(sub)

    TetraList=[]
    data=Dict["Tetra"]
    for cell, Typ in data:
        print cell
        aTetra=vtkTetra()
        for i in range(len(cell)):
            aTetra.GetPointIds().SetId(i,cell[i])
        aTetraGrid = vtk.vtkUnstructuredGrid()
        aTetraGrid.InsertNextCell(aTetra.GetCellType(), aTetra.GetPointIds())
        aTetraGrid.SetPoints(points)
        TetraList.append(aTetraGrid)

    return (SubLat, points, captions, color, TetraList)

def Plot(InputFile, HasCaption):
    # Read the data into a vtkPolyData using the functions in ReadPoints.py
    SubLat, points, captions, color, tetraList=Read(InputFile)

    data = vtk.vtkPolyData()
    data.SetPoints(points)
    data.GetPointData().SetScalars(color) 

    tActorList=[]
    for e in tetraList:
        tMap=vtkDataSetMapper()
        tMap.SetInputData(e)
        tActor=vtkActor()
        tActor.SetMapper(tMap)
        tActor.GetProperty().SetDiffuseColor(0, 1, 0)
        tActorList.append(tActor)

    ball = vtk.vtkSphereSource()
    ball.SetRadius(0.04)
    ball.SetThetaResolution(8)
    ball.SetPhiResolution(8)

    ballGlyph = vtk.vtkGlyph3D()
    if version<=5:
        ballGlyph.SetInput(data)
    else:
        ballGlyph.SetInputData(data)
    ballGlyph.SetSourceConnection(ball.GetOutputPort())
    ballGlyph.SetScaleModeToDataScalingOff()
    ballGlyph.SetColorModeToColorByScalar()
    ballGlyph.SetScaleFactor(1.0)


    colorTransferFunction = vtk.vtkColorTransferFunction()
    colorTransferFunction.AddRGBPoint(0, 1.0, 0.0, 0.0)
    colorTransferFunction.AddRGBPoint(1, 0.0, 0.0, 1.0)
    colorTransferFunction.AddRGBPoint(2, 0.0, 1.0, 0.0)
    colorTransferFunction.AddRGBPoint(3, 1.0, 1.0, 0.0)
    colorTransferFunction.AddRGBPoint(4, 0.0, 1.0, 1.0)
    colorTransferFunction.AddRGBPoint(5, 1.0, 1.0, 0.0)

    ballMapper = vtkPolyDataMapper()
    ballMapper.SetInputConnection(ballGlyph.GetOutputPort())
    ballMapper.SetLookupTable(colorTransferFunction)
    ballActor = vtkActor()
    ballActor.SetMapper(ballMapper)

    writer=vtkXMLPolyDataWriter()
    writer.SetFileName("test.vti")
    writer.SetInputConnection(ballGlyph.GetOutputPort())
    writer.Update()

    legend=vtkLegendBoxActor()
    legend.SetNumberOfEntries(4)
    for e in SubLat:
        legend.SetEntry(e, ball.GetOutput(), str(e),
                colorTransferFunction.GetColor(e)) 
    legend.BorderOff()
    legend.SetWidth(0.1)
    legend.SetHeight(0.1)
    legend.SetDisplayPosition(10,5)
    #txtProp=legend.GetEntryTextProperty()
# Create the Renderer, Window and Interator
    ren = vtkRenderer()
    ren.AddActor(ballActor)
    ren.AddActor(legend)
    for a in tActorList:
        ren.AddActor(a)
    ren.SetBackground(0.4, 0.4, 0.4)
    return ren

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="use file path to find the input file")
    parser.add_argument("-nc", "--nocaption", action="store_true", help="turn of caption")
    args = parser.parse_args()
    InputFile=os.path.abspath(args.file)

    ren=Plot(InputFile, not args.nocaption)
    renWin = vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetWindowName("Lattice")
    renWin.SetSize(1000, 800)
    iren = vtkRenderWindowInteractor()
    #iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
    iren.SetRenderWindow(renWin)

    iren.Initialize()
    iren.Start()
