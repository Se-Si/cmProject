import java.io.IOException;
import java.util.ArrayList;


public class NBody {

	//Initial Variable definitions
	// Gravitational Constant
	static double g = 1.0;

	//


	public static void main (String[] argv) throws IOException {

		// Number of timesteps
		int numstep = 5000;
		// Size of timestep
		double dt = 0.01;
		// Initial time
		double t = 0;

		//Read Planet info from SolarParticles.dat into 'particles' list

		ArrayList<Particle3D> particles = Particle3D.readFile("SolarParticles.dat");
		int particleCount = (particles.size() - 1);    	
		int tempParticleCount1;
		int tempParticleCount2;


		//VERLET START


		//////////////////////////

		//Force calculations, don't need to beware double-counting
		//NEED TO PREVENT CREATION OF SELF-ACTING FORCES
		//WILL REPLACE WITH FOR LOOPS
		//create gforcecurrent list and populate;
		//list structure:
		// 0-0, 0-1, ... 0-n, 1-0, 1-1, ..., n-n
		
		ArrayList<Vector3D> gforcecurrent = new ArrayList<>();
		
		tempParticleCount1 = 0;
		while (tempParticleCount1 <= particleCount){

			tempParticleCount2 = 0;
			while (tempParticleCount2 <= particleCount && te){

				gforcecurrent.add(gForce(particles.get(tempParticleCount1),particles.get(tempParticleCount2)));


				tempParticleCount2 += 1;

			}

			tempParticleCount1 += 1;

		}

		///////////////////////

		//Update Position of particles list based upon gforcecurrent list
		//sum gforcecurrent for each object, apply to object

		Vector3D tempGForce = new Vector3D (0.0,0.0,0.0);
		tempParticleCount1 = 0;
		while (tempParticleCount1 <= particleCount){

			tempParticleCount2 = 0;
			while (tempParticleCount2 <= particleCount){

				tempGForce = Vector3D.add(tempGForce,gforcecurrent.get(tempParticleCount2));
				tempParticleCount2 += 1;

			}

			particles.get(tempParticleCount1).leapPosition(tempGForce, dt);			
			tempParticleCount1 += 1;

		}

		////////////////////////
		

		//Force calculations, don't need to beware double-counting
		//NEED TO PREVENT CREATION OF SELF-ACTING FORCES
		//WILL REPLACE WITH FOR LOOPS
		//create gforcenew list and populate;
		//list structure:
		// 0-0, 0-1, ... 0-n, 1-0, 1-1, ..., n-n
		
		ArrayList<Vector3D> gforcenew = new ArrayList<>();
		
		tempParticleCount1 = 0;
		while (tempParticleCount1 <= particleCount){

			tempParticleCount2 = 0;
			while (tempParticleCount2 <= particleCount){

				gforcecurrent.add(gForce(particles.get(tempParticleCount1),particles.get(tempParticleCount2)));


				tempParticleCount2 += 1;

			}

			tempParticleCount1 += 1;

		}

		///////////////////////
		//test print statement, check if code runs
		System.out.printf("When I was a lad I served a term as office boy to an attorney's firm");

	}
	//END OF MAIN METHOD




	//System Static Methods

	//return GPE for two particles
	public static double potentialEnergy(Particle3D a, Particle3D b){

		Vector3D r = Particle3D.particleSeparation(a,b);
		double rmag = r.mag();
		double ma = a.getMass();
		double mb = b.getMass(); 	
		double pEnergy = g * ma * mb * (1/rmag);
		return pEnergy;  	

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
		gforce.multScalar(-1*ma*mb*g);
		gforce.divScalar(rmag*rmag);

		return gforce;

	}

}