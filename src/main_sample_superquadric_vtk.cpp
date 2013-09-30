#include <vtkVersion.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSuperquadricSource.h>
#include <vtkSuperquadric.h>
#include <vtkSampleFunction.h>
#include <vtkImageData.h>
#include <vtkContourFilter.h>

#include <pcl/io/pcd_io.h>

#include <pcl/console/print.h>
#include <pcl/console/parse.h>

using namespace pcl;

int
main (int argc,
      char **argv)
{
  PCL_INFO ("Syntax: %s [-e1 x.x] [-e2 y.y] [-a z.z] [-b t.t] [-c u.u] [-visualize 0/1] [-output xxx.pcd]\n", argv[0]);

  /// Parse command line arguments
  bool visualize_enabled = true;
  console::parse_argument (argc, argv, "-visualize", visualize_enabled);

  double epsilon_1 = 0.;
  console::parse_argument (argc, argv, "-e1", epsilon_1);

  double epsilon_2 = 0.;
  console::parse_argument (argc, argv, "-e2", epsilon_2);

  double a = 1.;
  console::parse_argument (argc, argv, "-a", a);

  double b = 1.;
  console::parse_argument (argc, argv, "-b", b);

  double c = 1.;
  console::parse_argument (argc, argv, "-c", c);

  std::string output_file = "";
  console::parse_argument (argc, argv, "-output", output_file);

  PCL_INFO ("### Superquadric parameters:\n   epsilon_1: %f\n   epsilon_2: %f\n   a: %f\n   b: %f\n   c: %f\n",
            epsilon_1, epsilon_2, a, b, c);


  /// Sample the superquadric to a pcd file
  if (output_file != "")
  {
    vtkSmartPointer<vtkSuperquadric> superquadric = vtkSmartPointer<vtkSuperquadric>::New ();
    superquadric->SetThetaRoundness (epsilon_1); /// epsilon_1
    superquadric->SetPhiRoundness (epsilon_2); /// epsilon_2
    superquadric->SetScale (a, b, c);
    superquadric->SetSize (1.);

    vtkSmartPointer<vtkSampleFunction> sample_function = vtkSmartPointer<vtkSampleFunction>::New ();
    sample_function->SetSampleDimensions (200, 200, 200);
    sample_function->SetImplicitFunction (superquadric);
    sample_function->SetModelBounds (-1, 1, -1, 1, -1, 1);


    vtkSmartPointer<vtkContourFilter> sample_contour = vtkSmartPointer<vtkContourFilter>::New ();
    sample_contour->SetInputConnection (sample_function->GetOutputPort ());
    sample_contour->SetValue (0, 0.0);
    sample_contour->Update ();
    vtkPolyData *poly_data = sample_contour->GetOutput ();

    PCL_INFO ("poly data cells %ld, points %ld\n",
              poly_data->GetNumberOfCells (), poly_data->GetNumberOfPoints ());


    /// Convert poly data to PCL format
    PointCloud<PointXYZ> cloud;
    cloud.resize (poly_data->GetNumberOfPoints ());
    double xyz[3];
    for (vtkIdType p_i = 0; p_i < poly_data->GetNumberOfPoints (); ++p_i)
    {
      poly_data->GetPoint (p_i, &xyz[0]);
      cloud[p_i].x = xyz[0];
      cloud[p_i].y = xyz[1];
      cloud[p_i].z = xyz[2];
    }
    io::savePCDFile (output_file, cloud, true);
  }


  /// Visualize the super quadric
  if (visualize_enabled)
  {
    // Create a superquadric
    vtkSmartPointer<vtkSuperquadricSource> superquadricSource = vtkSmartPointer<vtkSuperquadricSource>::New();
    superquadricSource->SetThetaRoundness (epsilon_2);
    superquadricSource->SetPhiRoundness (epsilon_1);
    superquadricSource->SetScale (a, b, c);
    superquadricSource->SetPhiResolution (50);
    superquadricSource->SetThetaResolution (50);
    superquadricSource->SetSize (1.);

    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> superquadricMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    superquadricMapper->SetInputConnection(superquadricSource->GetOutputPort());

    vtkSmartPointer<vtkActor> superquadricActor = vtkSmartPointer<vtkActor>::New();
    superquadricActor->SetMapper(superquadricMapper);

    vtkSmartPointer<vtkRenderer> sceneRenderer =  vtkSmartPointer<vtkRenderer>::New();

    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();

    renderWindow->AddRenderer(sceneRenderer);

    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    sceneRenderer->AddActor (superquadricActor);
    renderWindow->Render();
    renderWindowInteractor->Start();
  }

  return EXIT_SUCCESS;
}
