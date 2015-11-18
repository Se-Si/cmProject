import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Scanner;

/**
 * A class representing a point particle moving in three-dimensional space, with setters, getters and constructors.
 * Includes instance methods for time/position integration support and kinetic energy calculations,
 * and static methods for particle relative positioning and particle construction from file scanning.
 *
 * @author Sebastiaan van de Bund
 * @author James Maroulis
 * @author Sean Sirur
 *
 * @version "26/10/2015"
 */

public class Particle3D {
    /*
     * Particle properties
     */
    private String name;
    private double mass;
    private Vector3D position;
    private Vector3D velocity;
    
    /*
     * Setters and Getters
     */
    
    /** Get the name of a particle.
     *
     * @return a String representing the name of the particle.
     */
    public String getName() {
        return name;
    }
    
    /** Get the position of a particle.
     *
     * @return a Vector3D representing the position of the particle.
     */
    public Vector3D getPosition() { return position; }
    
    /** Get the velocity of a particle.
     *
     * @return a Vector3D representing the velocity of the particle.
     */
    public Vector3D getVelocity() { return velocity; }
    
    /** Get the mass of a particle.
     *
     * @return a double representing the mass.
     */
    public double getMass()     { return mass; }
	   
    /** Set the velocity of a particle.
     * @param v a Vector3D representing the particle velocity.
     */
    public void setVelocity(Vector3D v) { velocity = new Vector3D(v); }
    
    /** Set the position of a particle.
     *
     * @param p a Vector3D representing the particle position.
     */
    public void setPosition(Vector3D p) { position = new Vector3D(p); }
    
    /** Set the mass of a particle.
     *
     * @param m a double representing the mass of the particle.
     */
    public void setMass(double m)     { this.mass = m; }
    
    /** Set the name of a particle.
     *
     * @param name a String representing the particle name.
     */
    
    public void setName(String name) {
        this.name = name;
    }
    
    /*
     * Constructors
     */
    
    /** Default constructor. Sets all numerical particle properties to zero.
     *  Particle name will be set to null.
     */
    public Particle3D() {
        this.position = new Vector3D();
        this.velocity = new Vector3D();
        this.setName(null);
        this.setMass(0.0);
    }
    
    /** Explicit constructor. Constructs a new Particle1D with
     * explicitly given name, position, velocity, and mass.
     *
     * @param n a String that defines the name.
     * @param m a double that defines the mass.
     * @param p a Vector3D that defines the position.
     * @param v a Vector3D that defines the velocity.
     */
    public Particle3D(String n, double m, Vector3D p, Vector3D v) {
        this.position = new Vector3D();
        this.velocity = new Vector3D();
        this.setName(n);
        this.setMass(m);
        this.setPosition(p);
        this.setVelocity(v);
    }
    
    /** Returns a String representation of Particle3D, in VMD format as such:
     * label x y z
     *
     * @return the string representation.
     */
    public String toString() {
        return String.format("%s %f %f %f", name, getPosition().getX(), getPosition().getY(), getPosition().getZ());
    }
    
    /*
     * Instance Methods
     */
    
    /** Time integration support: evolve the velocity
     *  according to dv = f/m * dt.
     *
     * @param dt a double that is the timestep.
     * @param force a Vector3D that is the current force on the particle.
     */
    public void leapVelocity(Vector3D force, double dt){
        Vector3D f = new Vector3D(force);
        f.multScalar(dt);
        f.divScalar(this.mass);
        this.velocity = Vector3D.add(this.velocity, f);
    }
    
    /** Time integration support: evolve the position
     *  according to dx = v * dt.
     *
     * @param dt a double that is the timestep.
     */
    public void leapPosition(double dt){
        Vector3D v = new Vector3D(this.velocity);
        v.multScalar(dt);
        this.position = Vector3D.add(this.position, v);
    }
    
    /** Time integration support: evolve the position
     *  according to dx = v * dt + 0.5 * a * dt^2.
     *
     * @param dt a double that is the timestep.
     * @param force a Vector3D that is the current force.
     */
    public void leapPosition(Vector3D force, double dt){
        this.leapPosition(dt);
        Vector3D f = new Vector3D(force);
        f.divScalar(2.0*this.mass);
        f.multScalar(dt*dt);
        this.position = Vector3D.add(this.position, f);
    }
    
    /** Returns the kinetic energy of a Particle3D,
     *  calculated as 0.5*m*v^2.
     *
     * @return a double that is the kinetic energy.
     */
    public double kineticEnergy(){
        return this.getVelocity().magSquared() * this.getMass() * 0.5;
    }
    
    /*
     * Static Methods
     */
    
    /**
     * Reads the contents of a file and returns it in the form of an ArrayList of Particle3D objects.
     *
     * The file must contain 1 name and 7 numbers per line, and each line will give one particle:
     * The first number gives the mass.
     * The next three numbers give the initial position.
     * The last three numbers give the initial velocity.
     *
     * @param filename the name of the file from which the contents will be read
     */
    public static Particle3D[] readFile(String filename) throws FileNotFoundException {
        BufferedReader bufferedReader = new BufferedReader(new FileReader(filename));
        Scanner scanner = new Scanner(bufferedReader);

        //Read the first line, which is the number of particles
        int particleNum = Integer.parseInt(scanner.next());
        Particle3D[] particles = new Particle3D[particleNum];

        //Check whether there is a next line, read the name and 7 numbers, store them in a new Particle3D
        for(int i=0;i<particles.length;i++){
            particles[i] = new Particle3D(scanner.next(), scanner.nextDouble(), new Vector3D(scanner.nextDouble(), scanner.nextDouble(), scanner.nextDouble()), new Vector3D(scanner.nextDouble(), scanner.nextDouble(), scanner.nextDouble()));
        }
        scanner.close();
        return particles;
    }
    
    /**
     * Gives a Vector3D denoting the position of particle b relative to particle a.
     * @param a The first particle.
     * @param b The second particle.
     */
    public static Vector3D particleSeparation(Particle3D a, Particle3D b){
        return Vector3D.subtract(b.getPosition(),a.getPosition());
    }
}