import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Properties;
import java.util.Set;

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
 *
 * @version "24/11/2015"
 */
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
		Particle3D[] particles = Particle3D.readFile(argv[0]);
		//Array containing only the heliocentric bodies
		Particle3D[] heliocentricParticles = getHeliocentricBodies(particles);

		//Read the initial conditions from file
		Properties param = new Properties();
		param.load(new FileReader(argv[1]));
		Set<String> properties = param.stringPropertyNames();
		if (properties.contains("iterations")) {
			iterations = Integer.parseInt(param.getProperty("iterations"));
		}
		if (properties.contains("timestep")) {
			dt = Double.parseDouble(param.getProperty("timestep"));
		}
		if (properties.contains("gravconstant")) {
			g = Double.parseDouble(param.getProperty("gravconstant"));
		}

		//Create output file
		PrintWriter trajectoryOutput = new PrintWriter(new FileWriter(argv[2]));

		// Current forces at time t acting on all of the particles
		Vector3D[] currentForces;
		// New forces at time t+dt acting on all of the particles
		Vector3D[] newForces;

		/*
		 * Special indices and arrays for algorithms operating on specific bodies
		 */
		// Check for moon and assign index if present, otherwise set it to -1
		int lunaIndex = findParticle(particles, "Luna");
		// Check for moon and assign index if present
		int earthIndex = findParticle(particles, "Earth");
		// Store the index of the sun
		int solIndex = findParticle(particles, "Sol");

		/*
		 * Variables dealing with analysis of orbits, perihelions, aphelions and energy fluctuations
		 */
		// Store the indices of all heliocentric bodies
		double[] heliocentricOrbits = new double[heliocentricParticles.length];
		double[] lunarOrbit = new double[1];
		// Array containing the calculated values of the aphelions for all heliocentric particles
		double[] aphelions = new double[heliocentricParticles.length];
		// Array containing the calculated values of the perihelions for all heliocentric particles
		double[] perihelions = new double[heliocentricParticles.length];
		double lowestEnergy = Double.POSITIVE_INFINITY;
		double highestEnergy = Double.NEGATIVE_INFINITY;

		// Compute initial forces
		currentForces = totalInteractionForces(particles);

		//Write initial positions to file
		writePointsToFile(particles, 1, trajectoryOutput);


		/*
		 * Main algorithm
		 */
		for (int n=0; n<iterations; n++) {
			// Calculate and display progress
			if ((n+1) % (iterations / 20) == 0) {
				System.out.printf("Progress: %.0f%%\n", (float) n / (float) iterations * 100.f);
			}

			//Old and new data used for analysis of aphelions, perihelions and orbits
			double[] radialVelocityComponentsOld;
			double[] radialVelocityComponentsNew;
			Vector3D[] heliocentricPositionsOld;
			Vector3D[] heliocentricPositionsNew;
			Vector3D lunarPositionOld = new Vector3D();
			Vector3D lunarPositionNew = new Vector3D();

			// Compute radial velocity components before update
			radialVelocityComponentsOld = radialVelocityComponents(heliocentricParticles);
			// Compute positions before position update
			heliocentricPositionsOld = getPositions(heliocentricParticles);
			if(lunaIndex != -1) {
				lunarPositionOld = Vector3D.subtract(particles[lunaIndex].getPosition(), particles[earthIndex].getPosition());
			}

			// Leap all the particle positions
			leapPositions(particles, currentForces, dt);

			// Calculate new forces, with particles at time t+dt
			newForces = totalInteractionForces(particles);

			// Leap the velocities from the average of the forces
			leapVelocities(particles, elementAverage(currentForces, newForces), dt);

			// Forces at time t+dt become forces at time t for the next iteration
			currentForces = newForces;

			// Recompute the radial velocity components at time t+dt
			radialVelocityComponentsNew = radialVelocityComponents(heliocentricParticles);
			// Compute positions after position update, i.e. t+dt
			heliocentricPositionsNew = getPositions(heliocentricParticles);
			if(lunaIndex != -1) {
				lunarPositionNew = Vector3D.subtract(particles[lunaIndex].getPosition(), particles[earthIndex].getPosition());
			}

			incrementHeliocentricOrbits(heliocentricPositionsOld, heliocentricPositionsNew, heliocentricOrbits);
			if(lunaIndex != -1) {
				incrementLunarOrbit(lunarPositionOld, lunarPositionNew, lunarOrbit);
			}

			// Check whether apses has been passed
			checkApses(heliocentricParticles, radialVelocityComponentsOld, radialVelocityComponentsNew, aphelions, perihelions);

			// Calculate total energy and monitor energy fluctuations
			double totalEnergy = totalEnergy(particles);
			if (totalEnergy < lowestEnergy) {
				lowestEnergy = totalEnergy;
			}
			if (totalEnergy > highestEnergy) {
				highestEnergy = totalEnergy;
			}

			// Print output to file every 10 iterations
			if (n % 10 == 0) {
				writePointsToFile(particles, n + 1, trajectoryOutput);
			}

			// Increase time by timestep
			t += dt;
		}

		/*
		 * Display results of the analysis
		 */
		//Energy fluctuations, assuming the real energy is roughly the average of the lowest and highest energy
		System.out.println("\n--Energy Fluctuations--");
		System.out.printf("Global energy fluctuation: %.9f,  Relative error estimation: %f\n", Math.abs(highestEnergy-lowestEnergy),
				Math.abs(highestEnergy-lowestEnergy)/(Math.abs(highestEnergy+lowestEnergy)/2.0));

		// Aphelion and perihelion of heliocentric bodies
		System.out.println("\n--Aphelions and Perihelions--");
		for (int i=0; i<aphelions.length; i++) {
			System.out.printf("%s: ap = %f   pe = %f\n", heliocentricParticles[i].getName(), aphelions[i], perihelions[i]);
		}

		// Total partial orbits and orbital periods for heliocentric bodies and the moon
		System.out.println("\n--Number of Orbits and Periods--");
		for(int i=0; i<heliocentricOrbits.length;i++){
			System.out.printf("%s: %f orbits   -->   T = %f days = %f years\n", heliocentricParticles[i].getName(),
									heliocentricOrbits[i], t/heliocentricOrbits[i], t/heliocentricOrbits[i]/365.25);
		}
		if(lunaIndex != -1) {
			System.out.printf("%s: %f orbits   -->   T = %f days = %f years\n", particles[lunaIndex].getName(),
					lunarOrbit[0], t / lunarOrbit[0], t / lunarOrbit[0] / 365.25);
		}

		// Verification of Kepler's 3rd Law
		System.out.println("\n--Kepler's 3rd Law--");
		for(int i=0;i<heliocentricOrbits.length;i++){
			double a = (aphelions[i]+perihelions[i])/2.0;
			System.out.printf("%s: T^2/a^3 = %.3f days^2/AU^3\n", heliocentricParticles[i].getName(),
																t/heliocentricOrbits[i]*t/heliocentricOrbits[i]/(a*a*a));
		}

		// Print the ratio between the years of Mercury and Venus
		System.out.println("\n--Ratio of Mercury and Venus Periods--");
		System.out.printf("T_ven/T_mer = %f\n", heliocentricOrbits[2]/heliocentricOrbits[1]);


		trajectoryOutput.close();
	}


	/*
	 * Static methods
	 */
	/**
	 * Calculates the Gravitational Potential Energy due to the configuration of two Particle3D objects.
	 *
	 * @param a - Particle3D a
	 * @param b - Particle3D b
	 * @return Double representing GPE value.
	 */
	public static double potentialEnergy(Particle3D a, Particle3D b){

		Vector3D r = Particle3D.particleSeparation(a,b);
		double rmag = r.mag();
		double ma = a.getMass();
		double mb = b.getMass();
		return -g * ma * mb * (1.0/rmag);

	}

	/**
	 * Calculates the Newtonian Gravitational Force between two Particle3D objects.
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
	 * Calculates the Total Gravitational Force acting upon each body of an N-body Particle3D array due to all others.
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
	 * Calculates the total energy of an N-Body array of Particle3D objects as the sum of total kinetic and
	 * total potential energies.
	 * @param particles - Array of Particle3D objects.
	 * @return Double representing total particle energy.
	 */
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
	 * Convenience method; applies Particle3D.leapVelocity() to each object in a Particle3D array.
	 * Particles array element i will be acted upon by forces array element i.
	 *
	 * @param particles - Particle3D array.
	 * @param forces - Vector3D array of forces.
	 * @param dt - Timestep Value
	 */
	public static void leapVelocities(Particle3D[] particles, Vector3D[] forces, double dt){
		for(int i=0;i<particles.length;i++){
			particles[i].leapVelocity(forces[i], dt);
		}
	}

	/**
	 * Convenience method; applies Particle3D.leapPosition() to each object in a Particle3D array.
	 * particles array element i will be acted upon by forces array element i.
	 *
	 * @param particles - Particle3D array.
	 * @param forces - Vector3D array of forces.
	 * @param dt - Timestep Value.
	 */
	public static void leapPositions(Particle3D[] particles, Vector3D[] forces, double dt){
		for(int i=0;i<particles.length;i++){
			particles[i].leapPosition(forces[i], dt);
		}
	}

	/**
	 * Averages i elements of two equally sized Vector3D arrays, creating a new Vector3D array from the results.
	 * @param a - Vector3D array a.
	 * @param b - Vector3D array b.
	 * @return - Average Vector3D array.
	 */
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
		for(int i=0;i<radialVelocityComponentOld.length;i++){
			if((radialVelocityComponentOld[i] > 0.0) && (radialVelocityComponentNew[i] < 0.0)){
				aphelions[i] = particles[i].getPosition().mag();
			}
			if((radialVelocityComponentOld[i] < 0.0) && (radialVelocityComponentNew[i] > 0.0)){
				perihelions[i] = particles[i].getPosition().mag();
			}
		}
	}

	//Increments the partial orbits for heliocentric orbits
	public static void incrementHeliocentricOrbits(Vector3D[] oldPositions, Vector3D[] newPositions, double[] heliocentricOrbits){
		for(int i=0;i<heliocentricOrbits.length;i++){
			heliocentricOrbits[i] += Math.acos(Vector3D.dot(oldPositions[i], newPositions[i])
					/(oldPositions[i].mag()*newPositions[i].mag()))/(2.0*Math.PI);
		}
	}

	//Increments the partial orbits for the moon's geocentric orbit
	//lunarOrbit is only an array in order to allow it to be passed by reference
	public static void incrementLunarOrbit(Vector3D oldPosition, Vector3D newPosition, double[] lunarOrbit){
		lunarOrbit[0] += Math.acos(Vector3D.dot(oldPosition, newPosition)
				/(oldPosition.mag()*newPosition.mag()))/(2.0*Math.PI);
	}

    // Find the index of a planet or moon in the particles array, and if not found return -1
	public static int findParticle(Particle3D[] particles, String name) {
		for(int i=0;i<particles.length;i++){
			if(particles[i].getName().equals(name)){
				return i;
			}
		}
		return -1;
	}

	//Return an array of positions from an array of particles
	public static Vector3D[] getPositions(Particle3D[] particles){
		Vector3D[] positions = new Vector3D[particles.length];
		for(int i=0;i<particles.length;i++){
			positions[i] = particles[i].getPosition();
		}
		return positions;
	}

	// Return an array of all heliocentric bodies
	public static Particle3D[] getHeliocentricBodies(Particle3D[] particles){
		int heliocentricBodiesNum = 0;
		boolean[] isHeliocentrics = new boolean[particles.length];
		for(int i=0;i<particles.length;i++){
			if(isHeliocentric(particles[i].getName())){
				heliocentricBodiesNum++;
				isHeliocentrics[i] = true;
			} else {
				isHeliocentrics[i] = false;
			}
		}
		Particle3D[] heliocentricBodies = new Particle3D[heliocentricBodiesNum];
		int k=-1;
		for(int i=0;i<isHeliocentrics.length;i++){
			if(isHeliocentrics[i]){
				//Assign the reference so that their overlap refers to the same particles
				heliocentricBodies[++k] = particles[i];
			}
		}
		return heliocentricBodies;
	}

	//Helper method for getHeliocentricBodies()
	public static boolean isHeliocentric(String name){
		return name.equals("Mercury") || name.equals("Venus") || name.equals("Earth") || name.equals("Mars")
				|| name.equals("Jupiter") || name.equals("Saturn") || name.equals("Uranus") || name.equals("Neptune")
				|| name.equals("Pluto") || name.equals("1P/Halley");
	}
}