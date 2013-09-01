#include "segmentation_utils.h"

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>

#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_lib_io.h>


using namespace pcl;

void
sq::triangulizeContour (pcl::PointCloud<pcl::PointXYZ>::ConstPtr contour,
                        pcl::PolygonMesh &mesh)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New ();
  for (size_t i = 0; i < contour->size (); ++i)
    points->InsertNextPoint ((*contour)[i].x, (*contour)[i].y, (*contour)[i].z);

  vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New ();
  polygon->GetPointIds ()->SetNumberOfIds (contour->size ());
  for (size_t i = 0; i < contour->size (); ++i)
    polygon->GetPointIds ()->SetId (i, i);

  vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New ();
  polygons->InsertNextCell (polygon);

  vtkSmartPointer<vtkPolyData> polygon_poly_data = vtkSmartPointer<vtkPolyData>::New ();
  polygon_poly_data->SetPoints (points);
  polygon_poly_data->SetPolys (polygons);

  vtkSmartPointer<vtkTriangleFilter> triangle_filter = vtkSmartPointer<vtkTriangleFilter>::New ();
  triangle_filter->SetInput (polygon_poly_data);
  vtkSmartPointer<vtkPolyData> triangle_poly_data = vtkSmartPointer<vtkPolyData>::New ();
  triangle_filter->SetOutput (triangle_poly_data);
  triangle_filter->Update ();

  pcl::io::vtk2mesh (triangle_poly_data, mesh);
}
