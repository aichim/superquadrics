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

int main(int argc, char *argv[])
{
  // Create a superquadric
  vtkSmartPointer<vtkSuperquadricSource> superquadricSource = vtkSmartPointer<vtkSuperquadricSource>::New();
  superquadricSource->SetPhiRoundness(1.1);
  superquadricSource->SetThetaRoundness(.2);

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

  // Add actors to the renderers
  sceneRenderer->AddActor (superquadricActor);

  // Render once to figure out where the background camera will be
  renderWindow->Render();

  // Render again to set the correct view
  renderWindow->Render();

  // Interact with the window
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
