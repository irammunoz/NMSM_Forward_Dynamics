// Include OpenSim and functions
#include <OpenSim/OpenSim.h>
#include <OpenSim/Analyses/MuscleAnalysis.h>
#include <OpenSim/Analyses/ProbeReporter.h>
#include <OpenSim/Tools/AnalyzeTool.h>
#include "OpenSim/Common/STOFileAdapter.h"

// This allows us to use OpenSim functions, classes, etc., without having to
// prefix the names of those things with "OpenSim::".
using namespace OpenSim;

// This allows us to use SimTK functions, classes, etc., without having to
// prefix the names of those things with "SimTK::".
using namespace SimTK;

using namespace std;

class Arm22Torque : public OpenSim::Force {
    OpenSim_DECLARE_CONCRETE_OBJECT(Arm22Torque, Force);

public:
    /**
     * Constructor
     *
     * @param aModel Model to be controlled
     */
    Arm22Torque(const std::string& acoordName, double aq_upper, 
        double aq_lower, double aK1, double aK2, double aB) : 
        Force(), coordName(acoordName), q_upper(aq_upper), q_lower(aq_lower), K1(aK1), K2(aK2), B(aB)
    {
    }

    void computeForce(const SimTK::State& s,
        SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
        SimTK::Vector& generalizedForces) const override
    {
        SimTK::ReferencePtr<Coordinate> _coord;
        _coord = _model->updCoordinateSet().get(coordName);
        applyGeneralizedForce(s, *_coord, calcLimitForce(s,_coord), generalizedForces);

    }

    double Arm22Torque::calcLimitForce(const SimTK::State& s, SimTK::ReferencePtr<Coordinate> &_coord) const
    {

        double q = _coord->getValue(s);
        double qdot = _coord->getSpeedValue(s);

        double d = 0;
        double tau_p = 0;

        d = q - q_upper;

        // is angle above upper limit of ROM ?
        if (d > 0.0) {
            tau_p = -K2 * d * d;
        }

        // is angle below lower limit of ROM ?
        d = q - q_lower;
        if (d < 0.0) {
            tau_p = K2 * d * d;
        }

        // add a small amount of damping and overall stiffness
        tau_p = tau_p - B * qdot - K1 * q;

        return tau_p;

    }

    // This section contains the member variables of this controller class.
private:

    /**  */
    std::string coordName;
    double q_upper;
    double q_lower; 
    double K1; 
    double K2; 
    double B;
 
};

//______________________________________________________________________________
/**
 * This controller will try to track a desired trajectory of the block in
 * the Arm22 model.
 */
class Arm22Controller : public Controller {
    OpenSim_DECLARE_CONCRETE_OBJECT(Arm22Controller, Controller);

    // This section contains methods that can be called in this controller class.
public:
    /**
     * Constructor
     *
     * @param aModel Model to be controlled
     */
    Arm22Controller(void) : Controller()
    {
    }

    /**
     * This function is called at every time step for every actuator.
     *
     * @param s Current state of the system
     * @param controls Controls being calculated
     */
    void computeControls(const SimTK::State& s, SimTK::Vector& controls) const override
    {
        // Get the current time in the simulation.
        double t = s.getTime();

        Vector A(2); A[0] = 0.58; A[1] = 0.1;
        Vector tau(2); tau[0] = 1.5; tau[1] = 1.5;
        Vector tau_r(2); tau_r[0] = 0.5; tau_r[1] = 0.5;
        Vector T(2); T[0] = 3.0; T[1] = 3.0;
        Vector f0(2); f0[0] = 1 / T[0]; f0[1] = 1 / T[1];

        Vector t_offset(2); t_offset[0] = 0.0; t_offset[1] = 1.5;
        Vector u_offset(2); u_offset[0] = 0.32; u_offset[1] = 0.05;

        // Get pointers to each of the muscles in the model.
        auto TRImedMuscle = dynamic_cast<const Muscle*>  (&getActuatorSet().get(0));
        auto BRAMuscle = dynamic_cast<const Muscle*> (&getActuatorSet().get(1));

        double TRImedControl = TRImedMuscle->getMinControl(),
            BRAControl = BRAMuscle->getMinControl();

        Vector sumf(2);
        Vector u(2);

        sumf[0] = 0.0; sumf[1] = 0.0;
        u[0] = 0.0; u[1] = 0.0;

        Real d1 = 0;
        Real d2 = 0;
        Real d3 = 0;
        Real d4 = 0;
        Real d5 = 0;
        Real d6 = 0;

        Real dum1 = 0;
        Real dum2 = 0;


        for (int i = 0; i < 2; ++i) {
            for (int n = 1; n < 1000; ++n) {
                d1 = (sin(n * Pi * f0[i] * tau[i]) / (n * Pi * f0[i] * tau[i]));
                d2 = (sin(n * Pi * f0[i] * tau_r[i]) / (n * Pi * f0[i] * tau_r[i]));

                d3 = Pi * n * f0[i];
                d4 = (t + t_offset[i]);
                d5 = (tau[i] - tau_r[i]);

                d6 = cos(2 * d3 * d4 - d3 * d5);
                sumf[i] = sumf[i] + d1 * d2 * d6;
                dum1 = sumf[0];
                dum2 = sumf[1];
            }
            u[i] = u_offset[i] + A[i] * (tau[i] / T[i]) + 2 * A[i] * (tau[i] / T[i]) * sumf[i];
        }

        TRImedControl = u[1];
        BRAControl = u[0];

        // Don't allow any control value to be greater than one.
        if (TRImedControl > TRImedMuscle->getMaxControl())
            TRImedControl = TRImedMuscle->getMaxControl();
        if (BRAControl > BRAMuscle->getMaxControl())
            BRAControl = BRAMuscle->getMaxControl();

        // DeGrooteFregly muscle has only one control
        Vector muscleControl(1, TRImedControl);
        // Add in the controls computed for this muscle to the set of all model controls
        TRImedMuscle->addInControls(muscleControl, controls);
        // Specify control for other actuator (muscle) controlled by this controller
        muscleControl[0] = BRAControl;
        BRAMuscle->addInControls(muscleControl, controls);
    }

};

int main()
{
    try {

        // Create an OpenSim model from the model file provided.
        Model osimModel("arm22_dG.osim");
        osimModel.setUseVisualizer(true);
        
        // Create the controller.
        Arm22Controller* controller = new Arm22Controller();

        //Arm22Torque* force = new Arm22Torque();

        Arm22Torque* clfshoulder = new Arm22Torque("r_shoulder_elev",
            180 * Pi / 180, -90 * Pi / 180, 1, 10000, 1);
        osimModel.addForce(clfshoulder);

        Arm22Torque* clfelbow = new Arm22Torque("r_elbow_flex",
            130 * Pi / 180, 0 * Pi / 180, 1, 10000, 1);
        osimModel.addForce(clfelbow);

        // Give the controller the Model's actuators so it knows
        // to control those actuators.
        controller->setActuators(osimModel.updActuators());

        // Add the controller to the Model.
        osimModel.addController(controller);

        // Initialize the system and get the state representing the
        // system.
        State& si = osimModel.initSystem();
        // Fix the shoulder at its default angle and begin with the elbow flexed.
        osimModel.equilibrateMuscles(si);

        // Define non-zero (defaults are 0) states for the free joint.
        CoordinateSet& modelCoordinateSet =
            osimModel.updCoordinateSet();
        // Get the shoulder coordinate.
        Coordinate& shoulder = modelCoordinateSet.
            get("r_shoulder_elev");
        // Set shoulder position.
        shoulder.setValue(si, 0 * Pi / 180);
        // Get the elbow coordinate.
        Coordinate& elbow = modelCoordinateSet.
            get("r_elbow_flex");
        // Set elbow position.
        elbow.setValue(si, 0.349066);

        // Define the initial muscle states.
        const Set<Muscle>& muscleSet = osimModel.getMuscles();
        DeGrooteFregly2016Muscle* muscle1 = dynamic_cast<DeGrooteFregly2016Muscle*>(&muscleSet.get(0));
        DeGrooteFregly2016Muscle* muscle2 = dynamic_cast<DeGrooteFregly2016Muscle*>(&muscleSet.get(1));
        if ((muscle1 == NULL) || (muscle2 == NULL)) {
            throw OpenSim::Exception("ControllerExample: muscle1 or muscle2 is not an ActivationFiberLengthMuscle and example cannot proceed.");
        }
        muscle1->setActivation(si, 0.001); // muscle1 activation
        //muscle1->setFiberLength(si, 0.2); // muscle1 fiber length
        muscle2->setActivation(si, 0.5); // muscle2 activation
        //muscle2->setFiberLength(si, 0.2); // muscle2 fiber length

        // Configure the visualizer.
        osimModel.updMatterSubsystem().setShowDefaultGeometry(true);
        Visualizer& viz = osimModel.updVisualizer().updSimbodyVisualizer();
        viz.setBackgroundType(viz.SolidColor);
        viz.setBackgroundColor(Black);

        //ZYX(0,Pi/2,Pi/4) -> XYZ(-Pi/2, Pi/4, Pi/2)
        Vec3 camera_offset_direction(-Pi / 2, Pi / 4, Pi / 2);
        Vec3 camera_offset_distance(2, 2, 0); // meters
        Vec3 camera_look_at(0, 0, 0);

        SimTK::Transform camera_pose(camera_offset_distance);
        //camera_pose.updR().setRotationFromOneAxis(
        //	camera_offset_direction, YAxis);
        camera_pose.updR().setRotationToBodyFixedXYZ(camera_offset_direction);

        viz.setCameraTransform(camera_pose);
        //viz.pointCameraAt(camera_look_at, Vec3(0,1,0));

        // Setup ForceReporter and Manager
        ForceReporter* forces = new ForceReporter(&osimModel);
        osimModel.addAnalysis(forces);
        //osimModel.updAnalysisSet().adoptAndAppend(forces);

        // Create the manager for the simulation.
        Manager manager(osimModel);
        manager.setIntegratorMethod(Manager::IntegratorMethod::RungeKuttaMerson);
        //manager.setIntegratorAccuracy(1.0e-4);

        // Define the initial and final simulation times.
        double initialTime = 0.0;
        double finalTime = 5.0;

        double tsum = 0;
        SimTK::Vector time_lns(400);
        for (int i = 0; i < 400; ++i) {
            time_lns[i] = 0.0125;
            tsum = tsum + time_lns[i];
        }

        // Integrate from initial time to final time.
        si.setTime(initialTime);
        manager.initialize(si);
        
        manager.setIntegratorMinimumStepSize(0.010);
        manager.setIntegratorMaximumStepSize(0.010);
        std::cout << "\n\nIntegrating from " << initialTime
            << " to " << finalTime << std::endl;
        manager.integrate(finalTime);

        // Save the simulation results.
        auto controlsTable = osimModel.getControlsTable(); 
        STOFileAdapter_<double>::write(controlsTable, "Arm22_controls.sto");

        auto statesTable = manager.getStatesTable();
        STOFileAdapter_<double>::write(statesTable, "Arm22_states.sto");

        auto forcesTable = forces->getForcesTable();
        STOFileAdapter::write(forcesTable, "actuator_forces.sto");

        // Conduct an analysis using MuscleAnalysis and ProbeReporter.
        {
            // Create an AnalyzeTool setup file.
            AnalyzeTool analyze;
            analyze.setName("analyze");
            analyze.setModelFilename("arm22_dG.osim");
            analyze.setStatesFileName("Arm22_states.sto");
            analyze.updAnalysisSet().adoptAndAppend(new MuscleAnalysis());
            analyze.updAnalysisSet().adoptAndAppend(new ProbeReporter());
            analyze.setInitialTime(0.0);
            analyze.setFinalTime(5.0);
            analyze.updControllerSet().adoptAndAppend(
                new PrescribedController("Arm22_controls.sto"));
            analyze.print("Arm22_AnalyzeTool_setup.xml");
        }
        // Run the analysis
        AnalyzeTool analyze("Arm22_AnalyzeTool_setup.xml");
        analyze.run();

    }
    catch (const std::exception& ex) {

        // In case of an exception, print it out to the screen.
        std::cout << ex.what() << std::endl;

        // Return 1 instead of 0 to indicate that something
        // undesirable happened.
        return 1;
    }

    // If this program executed up to this line, return 0 to
    // indicate that the intended lines of code were executed.
    return 0;
}