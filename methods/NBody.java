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

		// Current forces at time t acting on all of the particles
		Vector3D[] currentForces;
		// New forces at time t+dt acting on all of the particles
		Vector3D[] newForces;

		// Variables used for analysis
		// Relative positions of planets from the sun
		Vector3D[] bodiesRelPost = new Vector3D[particles.length() - 1];
		// Number of orbts
		double[] numOrbits = new double[particles.length() - 1];
		// Array containing the calculated values of the aphelions for all the planets and pluto
		double[] aphelions = new double[particles.length];
		// Array containing the calculated values of the perihelions for all the planets and pluto
		double[] perihelions = new double[particles.length];
		double lowestEnergy = Double.POSITIVE_INFINITY;
		double highestEnergy = Double.NEGATIVE_INFINITY;
        // Index of the moon if included in the list of bodies
        static int moonIndex;
        // Index of the moon if included in the list of bodies
        static int earthIndex;

		// Compute initial forces
		currentForces = totalInteractionForces(particles);

		//Write initial positions to file
		writePointsToFile(particles, 1, trajectoryOutput);
		
		// Check for moon and assign index if present
		
		moonIndex = findMoon(particles);
        
        // Check for moon and assign index if present
		
		earthIndex = findEarth(particles);
        
        // Find start-of-system relative planet positions
        
        bodiesRelPos = setupPlanetsRelPos(particles, moonIndex);

        if (moonIndex != -1) {
        
            bodiesRelPos[moonIndex] = setupMoonRelPos(particles[moonIndex], particles[earthIndex]);
            
        }

        // Main Algorithm

		for(int n=0;n<iterations;n++) {
			//Calculate progress
			if(n%(iterations/20) == 0){
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

            // Update relative positions and number of orbits for each body
            numOrbits = planetOrbits(particles, bodiesRelPos, numOrbits, moonIndex);
            
            if (moonIndex != -1) {
            
                numOrbits[moonIndex] = moonOrbits(particles[moonIndex], particles[earthIndex], bodiesRelPos[moonIndex], numOrbits[moonIndex], moonIndex);
                
            }

			//TODO: add functionality to apply this to specific planets rather than every body in the system
			// Check whether apses has been passed
			checkApses(particles, radialVelocityComponentsOld, radialVelocityComponentsNew, aphelions, perihelions);

			// Print output to file every 10 iterations
			if(n%40 == 0) {
				writePointsToFile(particles, n + 1, trajectoryOutput);
			}

			// Increase time by timestep
			t += dt;
		}

		for(int i=0;i<aphelions.length;i++) {
			System.out.printf("r=%f\n",aphelions[i]);
		}
		
		for(int j=1;j<particles.length();j++) {
		    System.out.printf("Number of Orbits = %f\n", numOrbits[j]);

		System.out.println(Math.abs(highestEnergy-lowestEnergy));

		trajectoryOutput.close();
	}


	/*
	 * Static methods
	 */

	//return GPE for two particles
	public static double potentialEnergy(Particle3D a, Particle3D b){

		Vector3D r = Particle3D.particleSeparation(a,b);
		double rmag = r.mag();
		double ma = a.getMass();
		double mb = b.getMass();
		return -g * ma * mb * (1.0/rmag);

	}

	//Return gravitational force vector on a by b
	public static Vector3D gForce(Particle3D a, Particle3D b){
		Vector3D force;

		force = Particle3D.particleSeparation(a, b);
		double r = Particle3D.particleSeparation(a, b).mag();
		force.multScalar(g * a.getMass() * b.getMass() / (r * r * r));

		return force;
	}

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

	//Leap the velocities of an ArrayList of particles using an ArrayList of forces
	public static void leapVelocities(Particle3D[] particles, Vector3D[] forces, double dt){
		for(int i=0;i<particles.length;i++){
			particles[i].leapVelocity(forces[i], dt);
		}
	}

	//Leap the positions of an ArrayList of particles using an ArrayList of forces
	public static void leapPositions(Particle3D[] particles, Vector3D[] forces, double dt){
		for(int i=0;i<particles.length;i++){
			particles[i].leapPosition(forces[i], dt);
		}
	}

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

	//Write the coordinates of all particles to a file in VMD format
	public static void writePointsToFile(Particle3D[] particles, int pointNum, PrintWriter printWriter){
		printWriter.write(String.format("%d\n", particles.length));
		printWriter.write(String.format("Point = %d\n", pointNum));
		for(int i=0;i<particles.length;i++){
			printWriter.write(particles[i].toString());
		}
	}

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

	//Updates the array containing apsis values (passed by reference) by checking for a change of sign in the radial
	//velocity component
	public static void checkApses(Particle3D[] particles, double[] radialVelocityComponentOld, double[] radialVelocityComponentNew
														, double[] aphelions, double[] perihelions){
		boolean[] isAphelion = new boolean[radialVelocityComponentOld.length];
		boolean[] isPerihelion = new boolean[radialVelocityComponentOld.length];

		for(int i=0;i<radialVelocityComponentOld.length;i++){
			if((radialVelocityComponentOld[i] > 0.0) && (radialVelocityComponentNew[i] < 0.0)){
				aphelions[i] = particles[i].getPosition().mag();
			}
			if((radialVelocityComponentOld[i] < 0.0) && (radialVelocityComponentNew[i] > 0.0)){
				perihelions[i] = particles[i].getPosition().mag();
			}
		}
	}
	
	 /*

        Run this before the main algorithm starts
        It finds the relative positions of the planets at the start of the system
        
    */

    public static Vector3D setupPlanetsRelPos(Particle3D[] planets, int moonIndex) {

        private int loopOrbits;

        private Vector3D planetsRelPos = new Vector3D[planets.length()- 1];

        for (loopOrbits = 1; loopOrbits < planets[].length(); loopOrbits++) {

            if (loopOrbits != moonIndex) {
            
                planetsRelPos[i-1] = Vector3D.subtract(planets[i].getPosition(), planets[0].getPosition());
                
            }
            
        }
        
        return planetsRelPos;
        
    }
    
     /*

       This method assumes that the array planetsPreRelPos is set with all the starting positions of
       the planets and comets before the main algorithm's loop is started.
        
    */   
            
    public static double planetOrbits(Particle3D[] planets, Vector3D[] planetsPreRelPos, double[] numPlanetOrbits, int moonIndex) {

        private int loopOrbits;
        private Vector3D[] planetsCurrentRelPos = new Vector3D[planets.length() - 1];
        
        
        for (loopOrbits = 1; loopOrbits < planets.length; loopOrbits++) {
        
            if (loopOrbits != moonIndex) {
        
                planetsCurrentRelPos[i-1] = Vector3D.subtract(planets[i].getPosition(), planets[0].getPosition());
                
                numPlanetOrbits[i-1] = numPlanetOrbits[i-1] + java.lang.Math.acos(Vector3D.dot(planetsPreRelPos[i-1].getPosition(), planetsCurrentRelPos[i-1].getPosition()) / (planetsPreRelPos[i-1].mag * planetsCurrentRelPos[i-1].mag * 2 * Math.PI));
                
                planetsPreRelPos = new Vector3D(planetsCurrentRelPos);
         
            }
                
        }
        
        return numPlanetOrbits;
      
    }
    
    /*

        This is the pre-algorithm setup method for the relative starting position of the moon to the earth.
        
    */
        
    public static Vector3D setupMoonRelPos(Particle3D moon, Particle3D earth) {

        private Vector3D moonRelPos = new Vector3D();

        moonRelPos= Vector3D.subtract(moon.getPosition(), earth.getPosition());
        
        return moonRelPos;
        
    }

    /*

        This is the main algorithm method for calculating the number of orbits and updating the 
        relative position of the moon to the earth
        
    */

    public static double moonOrbits(Particle3D moon, Particle3D earth, Vector3D moonPreRelPos, double numMoonOrbits, int moonIndex) {

        private Vector3D moonCurrentRelPos = new Vector3D();

        moonCurrentRelPos = Vector3D.subtract(moon.getPosition(), earth.getPosition());
        
        numMoonOrbits = numMoonOrbits + java.lang.Math.acos(Vector3D.dot(moonPreRelPos.getPosition(), moonCurrentRelPos.getPosition()) / (moonPreRelPos.mag * moonCurrentRelPos.mag * 2 * Math.PI));
        
        moonPreRelPos = new Vector3D(moonCurrentRelPos);
            
        return numMoonOrbits;
            
    }    
    
    /*

        Use this to check for the moon and assign its index.

    */

    public static int findMoon(Particle3D[] particles) {

        for (loopOrbits = 0; loopOrbits < planets.length; loopOrbits++) {

            if (moon.getName() == "Luna") {

                moonIndex = loopOrbits;
                
            } else {
            
                moonIndex = -1;
                
            }
            
            return moonIndex;
            
        }
    
    }
    
    /*

        Use this to check for the earth and assign its index.

    */
    
    public static int findEarth(Particle3D[] particles) {

        for (loopOrbits = 0; loopOrbits < planets.length; loopOrbits++) {

            if (moon.getName() == "Earth") {

                earthIndex = loopOrbits;
                
            } else {
            
                earthIndex = -1;
                
            }
            
            return earthIndex;
            
        }
    
    }
        
}
