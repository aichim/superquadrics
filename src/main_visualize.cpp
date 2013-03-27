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

#include <pcl/console/print.h>
#include <pcl/console/parse.h>

using namespace pcl;

int
main (int argc,
      char **argv)
{
  PCL_INFO ("Syntax: %s [-e1 x.x] [-e2 y.y] [-a z.z] [-b t.t] [-c u.u]\n", argv[0]);

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

  PCL_INFO ("### Superquadric parameters:\n   epsilon_1: %f\n   epsilon_2: %f\n   a: %f\n   b: %f\n   c: %f\n",
            epsilon_1, epsilon_2, a, b, c);



  if (visualize_enabled)
  {
    // Create a superquadric
    vtkSmartPointer<vtkSuperquadricSource> superquadricSource = vtkSmartPointer<vtkSuperquadricSource>::New();
    /// Piet Hein Super-Egg
    //  superquadricSource->SetThetaRoundness (1.);
    //  superquadricSource->SetPhiRoundness (0.8);
    //  superquadricSource->SetScale (3., 3., 4.);

    superquadricSource->SetThetaRoundness (epsilon_1); /// epsilon_1
    superquadricSource->SetPhiRoundness (epsilon_2); /// epsilon_2
    superquadricSource->SetScale (a, b, c);

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
