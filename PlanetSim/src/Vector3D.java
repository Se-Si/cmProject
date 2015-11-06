/**
 * @author Sebastiaan van de Bund
 * @author James Maroulis
 * @author Sean Sirur
 */

public class Vector3D{
    private double x,y,z;
    
    /*
     * Constructors
     */
    /**
     * Default constructor. Creates a vector with all 3 components set to 0.
     */
    public Vector3D(){
        this.x=0;
        this.y=0;
        this.z=0;
    }
    /**
     * Explicit constructor. Creates a vector with all 3 components explicitly given.
     * @param xx a double as the x component
     * @param yy a double as the y component
     * @param zz a double as the z component
     */
    public Vector3D(double xx, double yy, double zz){
        this.x = xx;
        this.y = yy;
        this.z = zz;
    }
    /**
     * Copy constructor. Constructs a new Vector3D by copying the components of another Vector3D instance.
     * @param v the original Vector3D to be copied
     */
    public Vector3D(Vector3D v){
        this.x = v.getX();
        this.y = v.getY();
        this.z = v.getZ();
    }


     /*
     * Setters and getters
     */
    /**
     * Set the x component of the vector.
     * @param xx a double as the x component
     */
    public void setX(double xx){ this.x = xx; }
    /**
     * Set the y component of the vector.
     * @param yy a double as the y component
     */
    public void setY(double yy){ this.y = yy; }
    /**
     * Set the z component of the vector.
     * @param zz a double as the z component
     */
    public void setZ(double zz){ this.z = zz; }
    /**
     * Gets the x component of the vector.
     */
    public double getX(){ return this.x; }
    /**
     * Gets the y component of the vector.
     */
    public double getY(){ return this.y; }
    /**
     * Gets the z component of the vector.
     */
    public double getZ(){ return this.z; }


    /*
     * Instance methods
     */
    /**
     * Gives the square of the magnitude of the vector.
     */
    public double magSquared(){
        return x*x+y*y+z*z;
    }
    /**
     * Gives the magnitude of the vector.
     */
    public double mag(){
        return Math.sqrt(magSquared());
    }
    /**
     * Multiplies the vector by a scalar.
     * @param s a scalar with which the vector is multiplied
     */
    public void multScalar(double s){
        this.x*=s;
        this.y*=s;
        this.z*=s;
    }
    /**
     * Divides the vector by a scalar.
     * @param s a scalar with which the vector is divided
     */
    public void divScalar(double s){
        multScalar(1.f/s);
    }
    /**
     * Gives a String representation of the vector.
     */
    public String toString(){
        String s = String.format("(%f, %f, %f)",x,y,z);
        return s;
    }

    /*
     * Static methods
     */
    /**
     * Adds 2 vectors together.
     * @param v1 the first vector to be added
     * @param v2 the second vector to be added
     */
    public static Vector3D add(Vector3D v1, Vector3D v2){
        return new Vector3D(v1.getX()+v2.getX(),v1.getY()+v2.getY(),v1.getZ()+v2.getZ());
    }
    /**
     * Subtracts a vector v2 from a vector v1
     * @param v1 the first vector
     * @param v2 the second vector
     */
    public static Vector3D subtract(Vector3D v1, Vector3D v2){
        return new Vector3D(v1.getX()-v2.getX(),v1.getY()-v2.getY(),v1.getZ()-v2.getZ());
    }
    /**
     * Calculates the dot product of vectors v1 and v2
     * @param v1 the first vector
     * @param v2 the second vector
     */
    public static double dot(Vector3D v1, Vector3D v2){
        return v1.getX()*v2.getX()+v1.getY()*v2.getY()+v1.getZ()*v2.getZ();
    }
    /**
     * Calculates the cross product of vectors v1 and v2
     * @param v1 the first vector
     * @param v2 the second vector
     */
    public static Vector3D cross(Vector3D v1, Vector3D v2){
        return new Vector3D(v1.getY()*v2.getZ() - v1.getZ()*v2.getY(),
                            v1.getZ()*v2.getX() - v1.getX()*v2.getZ(),
                            v1.getX()*v2.getY() - v1.getY()*v2.getX());
    }
    /**
     * Checks whether the vectors v1 and v2 are equivalent by comparing components.
     * @param v1 the first vector
     * @param v2 the second vector
     */
    public static boolean isEquivalent(Vector3D v1, Vector3D v2){
        //Threshold value
        final double EPSILON = 1E-10;
        return (Math.abs(v1.getX() - v2.getX()) < EPSILON) && (Math.abs(v1.getY() - v2.getY()) < EPSILON) && (Math.abs(v1.getZ() - v2.getZ()) < EPSILON);
    }
}
