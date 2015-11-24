/** 
 * This code is to simulate the orbital motion of an N-body system due to gravitational forces,
 * with functionality in place to track numbers of orbits, orbital period,
 * and apoapsis/periapsis of each body, each assumed to be exhibiting orbital behaviour.
 * 
 * Included are several methods for the calculation of relevant physical quantities,
 * and several convenience methods for the application of various Vector3D and Particle3D methods to array formats.
 * 
 * @author Sebastiaan van de Bund
 * @author James Maroulis
 * @author Sean Sirur
 */
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Properties;
import java.util.Set;


public class NBody {
	// Simulation parameters, defaulted to 0
	// Number of timesteps
	static int iterations = 0;
	// Size of timestep
	static double dt = 0.0;
	// Gravitational Constant
	static double g = 0.0;

	public static void main (String[] argv) throws IOException {
		// Initial time
		double t = 0.0;


		//Read the particles and initial conditions from file
		//TODO: figure out a nice format for particles and maybe change the particle3D method a bit
		Particle3D[] particles = Particle3D.readFile(argv[0]);

		//Read the initial conditions from file
		Properties param = new Properties();
		param.load(new FileReader(argv[1]));
		Set<String> properties = param.stringPropertyNames();
		if(properties.contains("iterations")){
			iterations = Integer.parseInt(param.getProperty("iterations"));
		}
		if(properties.contains("timestep")){
			dt = Double.parseDouble(param.getProperty("timestep"));
		}
		if(properties.contains("gravconstant")){
			g = Double.parseDouble(param.getProperty("gravconstant"));
		}

		//Create output file
		PrintWriter trajectoryOutput = new PrintWriter(new FileWriter(argv[2]));

		// Variables used for analysis
		// Array containing the calculated values of the aphelions for all the planets and pluto
		double[] aphelions = new double[particles.length];
		for(int i=0;i<aphelions.length;i++){
			aphelions[i]=0.0;
		}
		// Array containing the calculated values of the perihelions for all the planets and pluto
		double[] perihelions = new double[particles.length];
		for(int i=0;i<aphelions.length;i++){
			perihelions[i]=0.0;
		}
		double lowestEnergy = Double.POSITIVE_INFINITY;
		double highestEnergy = Double.NEGATIVE_INFINITY;

		// Current forces at time t acting on all of the particles
		Vector3D[] currentForces;
		// New forces at time t+dt acting on all of the particles
		Vector3D[] newForces;

		// Compute initial forces
		currentForces = totalInteractionForces(particles);

		//Write initial positions to file
		writePointsToFile(particles, 1, trajectoryOutput);

		for(int n=0;n<iterations;n++) {
			//Calculate progress
			if(n%(iterations/10) == 0){
				System.out.printf("Progress: %.0f%%\n", (float)n/(float)iterations * 100.f);
			}

			double[] radialVelocityComponentsOld;
			double[] radialVelocityComponentsNew;

			// Compute radial velocity components before update
			radialVelocityComponentsOld = radialVelocityComponents(particles);

			// Leap all the particle positions
			leapPositions(particles, currentForces, dt);

			// Calculate new forces, with particles at time t+dt
			newForces = totalInteractionForces(particles);

			// Leap the velocities from the average of the forces
			leapVelocities(particles, elementAverage(currentForces, newForces), dt);

			// Forces at time t+dt become forces at time t for the next iteration
			currentForces = newForces;

			// Recompute the radial velocity components at time t+dt
			radialVelocityComponentsNew = radialVelocityComponents(particles);

			// Calculate total energy and monitor energy fluctuations
			double totalEnergy = totalEnergy(particles);
			if(totalEnergy < lowestEnergy){
				lowestEnergy = totalEnergy;
			}
			if(totalEnergy > highestEnergy){
				highestEnergy = totalEnergy;
			}

			//TODO: add functionality to apply this to specific planets rather than every body in the system
			// Check whether apses has been passed
			checkApses(particles, radialVelocityComponentsOld, radialVelocityComponentsNew, aphelions, perihelions);

			// Print output to file every 10 iterations
			if(n%20 == 0) {
				writePointsToFile(particles, n + 1, trajectoryOutput);
			}

			// Increase time by timestep
			t += dt;
		}

		for(int i=0;i<aphelions.length;i++) {
			System.out.printf("aphelion = %f,   perihelion = %f\n", aphelions[i], perihelions[i]);
		}

		System.out.println(Math.abs(highestEnergy-lowestEnergy)/Math.abs((highestEnergy+lowestEnergy)/2.0));

		trajectoryOutput.close();
	}


	/*
	 * Static methods
	 */
	 
	/**
	 * Calculates Gravitational Potential Energy of two Particle3D objects.
	 * 
	 * @param a - Particle3D a
	 * @param b - Particle3D b
	 * @return Double representing GPE value.
	 */
	//return GPE for two particles
	public static double potentialEnergy(Particle3D a, Particle3D b){

		Vector3D r = Particle3D.particleSeparation(a,b);
		double rmag = r.mag();
		double ma = a.getMass();
		double mb = b.getMass();
		return -g * ma * mb * (1.0/rmag);

	}
	
	/**
	 * Calculates Gravitational Force between two Particle3D objects.
	 * 
	 * @param a - Particle3D a
	 * @param b - Particle3D b
	 * @return Vector3D representing g-force of b upon a.
	 */
	//Return gravitational force vector on a by b
	public static Vector3D gForce(Particle3D a, Particle3D b){
		Vector3D force;

		force = Particle3D.particleSeparation(a, b);
		double r = Particle3D.particleSeparation(a, b).mag();
		force.multScalar(g * a.getMass() * b.getMass() / (r * r * r));

		return force;
	}

	/**
	 * Calculates the Total Gravitational Force acting upon each body of an N-body PArticle3D array due to those N bodies.
	 * 
	 * @param particles - Array of Particle3D objects.
	 * @return Vector3D array; Array element i represents gforce vector upon element i in particles array.
	 */
	public static Vector3D[] totalInteractionForces(Particle3D[] particles){
		Vector3D[] forces = new Vector3D[particles.length];

		// Calculate the force on the i-th particle
		for(int i=0;i<particles.length;i++){
			// Total force on i-th particle
			Vector3D force_i = new Vector3D(0.0, 0.0, 0.0);
			// Add forces due to all other particles except itself
			for(int j=0;j<particles.length;j++){
				if(i != j) {
					force_i = Vector3D.add(force_i, gForce(particles[i], particles[j]));
				}
			}
			forces[i] = force_i;
		}

		return forces;
	}

	/**
	 * Calculates the total energy of an N-Body array of Particle3D objects.
	 * @param particles - Array of Particle3D objects.
	 * @return Double representing total particle energy.
	 */
	//Return the total energy of the system of particles
	public static double totalEnergy(Particle3D[] particles){
		double totalKinetic = 0.0;
		double totalPotential = 0.0;

		//Calculate total kinetic energy
		for(int i=0;i<particles.length;i++){
			totalKinetic += particles[i].kineticEnergy();
		}

		//Calculate total potential energy
		for(int i=0;i<particles.length;i++){
			for(int j=0;j<i-1;j++){
				totalPotential+= -1.0 * g * particles[i].getMass() * particles[j].getMass()
										/ Particle3D.particleSeparation(particles[i], particles[j]).mag();
			}
		}

		return totalKinetic + totalPotential;
	}

	/**
	 * Convenience method; applies Particle3D.leapVelocity to each object in a Particle3D array.
	 * particles array element i will be acted upon by forces array element i.
	 * 
	 * @param particles - Particle3D array.
	 * @param forces - Vector3D array of forces.
	 * @param dt - Timestep Value
	 */
	//Leap the velocities of an ArrayList of particles using an ArrayList of forces
	public static void leapVelocities(Particle3D[] particles, Vector3D[] forces, double dt){
		for(int i=0;i<particles.length;i++){
			particles[i].leapVelocity(forces[i], dt);
		}
	}

	/**
	 * Convenience method; applies Particle3D.leapPosition to each object in a Particle3D array.
	 * particles array element i will be acted upon by forces array element i.
	 * 
	 * @param particles - Particle3D array.
	 * @param forces - Vector3D array of forces.
	 * @param dt - Timestep Value.
	 */
	//Leap the positions of an ArrayList of particles using an ArrayList of forces
	public static void leapPositions(Particle3D[] particles, Vector3D[] forces, double dt){
		for(int i=0;i<particles.length;i++){
			particles[i].leapPosition(forces[i], dt);
		}
	}

	/**
	 * Averages i elements of two Vector3D arrays, creates a new Vector3D array from the results.
	 * @param a - Vector3D array a.
	 * @param b - Vector3D array b.
	 * @return - Average Vector3D array.
	 */
	//Convenience method for computing the element-wise average of two equally sized ArrayLists
	public static Vector3D[] elementAverage(Vector3D[] a, Vector3D[] b){
		if(a.length == b.length) {
			Vector3D[] averages = new Vector3D[a.length];
			for (int i=0; i<a.length;i++){
				Vector3D vec = Vector3D.add(a[i], b[i]);
				vec.divScalar(2.0);
				averages[i] = vec;
			}
			return averages;
		} else {
			return new Vector3D[0];
		}
	}

//TODO: Finish JavaDoc from here onwards
	/**
	 * 
	 * @param particles
	 * @param pointNum
	 * @param printWriter
	 */
	//Write the coordinates of all particles to a file in VMD format
	public static void writePointsToFile(Particle3D[] particles, int pointNum, PrintWriter printWriter){
		printWriter.write(String.format("%d\n", particles.length));
		printWriter.write(String.format("Point = %d\n", pointNum));
		for(int i=0;i<particles.length;i++){
			printWriter.write(particles[i].toString());
		}
	}

	/**
	 * 
	 * @param particles
	 * @return
	 */
	//Computes the radial velocity components compared to the origin for an array of particles
	public static double[] radialVelocityComponents(Particle3D[] particles){
		Vector3D radialVector;
		double[] radialVelocityComponents = new double[particles.length];
		for(int i=0;i<particles.length;i++) {
			radialVector = new Vector3D(particles[i].getPosition());
			radialVector.divScalar(particles[i].getPosition().mag());
			radialVelocityComponents[i] = Vector3D.dot(particles[i].getVelocity(), radialVector);
		}
		return radialVelocityComponents;
	}

	/**
	 * 
	 * @param particles
	 * @param radialVelocityComponentOld
	 * @param radialVelocityComponentNew
	 * @param aphelions
	 * @param perihelions
	 */
	//Updates the array containing apsis values (passed by reference) by checking for a change of sign in the radial
	//velocity component
	public static void checkApses(Particle3D[] particles, double[] radialVelocityComponentOld, double[] radialVelocityComponentNew
														, double[] aphelions, double[] perihelions){
		boolean[] isAphelion = new boolean[radialVelocityComponentOld.length];
		boolean[] isPerihelion = new boolean[radialVelocityComponentOld.length];

		for(int i=0;i<radialVelocityComponentOld.length;i++){
			if((radialVelocityComponentOld[i] > 0.0) && (radialVelocityComponentNew[i] < 0.0)){
				if(aphelions[i] == 0) {
					aphelions[i] = particles[i].getPosition().mag();
				}
			}
			if((radialVelocityComponentOld[i] < 0.0) && (radialVelocityComponentNew[i] > 0.0)){
				if(perihelions[i] == 0) {
					perihelions[i] = particles[i].getPosition().mag();
				}
			}
		}
	}
}
