import java.io.IOException;
import java.util.ArrayList;


public class NBody {

	//Initial Variable definitions
	// Gravitational Constant
	static double g = 1.0;

	public static void main (String[] argv) throws IOException {

		// Number of timesteps
		int numstep = 5000;
		// Size of timestep
		double dt = 0.01;
		// Initial time
		double t = 0;

		//Read Planet info from SolarParticles.dat into 'particles' list
		ArrayList<Particle3D> particles = Particle3D.readFile("SolarParticles.dat");


		// Current forces at time t acting on all of the particles
		ArrayList<Vector3D> currentForces;
		// New forces at time t+dt acting on all of the particles
		ArrayList<Vector3D> newForces;


		//VERLET START

		// Compute initial forces
		currentForces = totalInteractionForces(particles);

		for(int n=0;n<numstep;n++) {
			// Leap all the particle positions
			leapPositions(particles, currentForces, dt);

			// Calculate new forces, with particles at time t+dt
			newForces = totalInteractionForces(particles);

			// Leap the velocities from the average of the forces
			leapVelocities(particles, elementAverage(currentForces, newForces), dt);
		}

	}


	//System Static Methods

	//return GPE for two particles
	public static double potentialEnergy(Particle3D a, Particle3D b){

		Vector3D r = Particle3D.particleSeparation(a,b);
		double rmag = r.mag();
		double ma = a.getMass();
		double mb = b.getMass();
		return -g * ma * mb * (1.0/rmag);

	}

	//return gravitational force vector of a on b
	public static Vector3D gForce(Particle3D a, Particle3D b){

		//Create relative position vector and relative direction vector
		Vector3D r = Particle3D.particleSeparation(a,b);
		double rmag = r.mag();
		Vector3D rhat = r;
		rhat.divScalar(rmag);

		//Create gravitational force for particle a acting on b
		Vector3D gforce = rhat;
		double ma = a.getMass();
		double mb = b.getMass();
		gforce.multScalar(-1.0*ma*mb*g);
		gforce.divScalar(rmag*rmag);

		return gforce;

	}

	public static ArrayList<Vector3D> totalInteractionForces(ArrayList<Particle3D> particles){
		ArrayList<Vector3D> forces = new ArrayList<>(particles.size());

		// Calculate the force on the i-th particle
		for(int i=0;i<particles.size();i++){
			// Total force on i-th particle
			Vector3D force_i = new Vector3D(0.0, 0.0, 0.0);
			// Add forces due to all other particles except itself
			for(int j=0;j<particles.size();j++){
				if(i != j) {
					Vector3D.add(force_i, gForce(particles.get(i), particles.get(j)));
				}
			}
			forces.set(i, force_i);
		}

		return forces;
	}

	//Return the total energy of the system of particles
	public static double totalEnergy(ArrayList<Particle3D> particles){
		double totalKinetic = 0.0;
		double totalPotential = 0.0;

		//Calculate total kinetic energy
		for(int i=0;i<particles.size();i++){
			totalKinetic += particles.get(i).kineticEnergy();
		}

		//Calculate total potential energy
		for(int i=0;i<particles.size();i++){
			for(int j=0;j<i-1;j++){
				totalPotential+= -1.0 * g * particles.get(i).getMass() * particles.get(j).getMass()
										/ Particle3D.particleSeparation(particles.get(i), particles.get(j)).mag();
			}
		}

		return totalKinetic + totalPotential;
	}

	//Leap the velocities of an ArrayList of particles using an ArrayList of forces
	public static void leapVelocities(ArrayList<Particle3D> particles, ArrayList<Vector3D> forces, double dt){
		for(int i=0;i<particles.size();i++){
			particles.get(i).leapVelocity(forces.get(i), dt);
		}
	}

	//Leap the positions of an ArrayList of particles using an ArrayList of forces
	public static void leapPositions(ArrayList<Particle3D> particles, ArrayList<Vector3D> forces, double dt){
		for(int i=0;i<particles.size();i++){
			particles.get(i).leapPosition(forces.get(i), dt);
		}
	}

	//Convenience method for computing the element-wise average of two equally sized ArrayLists
	public static ArrayList<Vector3D> elementAverage(ArrayList<Vector3D> a, ArrayList<Vector3D> b){
		if(a.size() == b.size()) {
			ArrayList<Vector3D> averages = new ArrayList<>();
			for (int i=0; i<a.size();i++){
				Vector3D vec = Vector3D.add(a.get(i), b.get(i));
				vec.divScalar(2.0);
				averages.set(i, vec);
			}
			return averages;
		} else {
			return new ArrayList<>();
		}
	}
}