using System.Collections.Generic;
using DelaunatorSharp;
using UnityEditor;
using UnityEngine;

public class SARSV2
{
    private static AudioReverbZone[] m_ReverbZones;
    private static GameObject m_CurrentCentrePoint;
    private const float SCALEFACTOR = 0.1f;

    [MenuItem("Tools/SARS V2.0")]
    private static void main()
    {
        m_ReverbZones = Object.FindObjectsOfType<AudioReverbZone>();

        foreach (var reverbZone in m_ReverbZones)
        {
            float totalSurfaceArea;
            float totalVolume;
            float[] avAbsorptionSpec;
            float avSmoothness;

            getCentrePoint(reverbZone);
            List<Vector3> points = getPoints(65535); // num points|| max num = 65535 due to int variable type

            int[] allTriangles = DelaunayTriangulation(points).ToArray();

            IList<RaycastHit> raycastCollisions = raycasts(points);      
            List<Vector3> CollisionPoints = rcCollisionPoints(raycastCollisions);

            // Debug Mesh
            Mesh mesh = new Mesh();
            mesh.SetVertices(CollisionPoints);
            mesh.SetTriangles(allTriangles, 0);
            m_CurrentCentrePoint.GetComponent<MeshFilter>().mesh = mesh;

            List<int> listAllTriangles = new List<int>(allTriangles);
            geometricCalculations(listAllTriangles, CollisionPoints, raycastCollisions, out avAbsorptionSpec, out totalVolume, out totalSurfaceArea, out avSmoothness);
            setParameters(reverbZone, avAbsorptionSpec, totalVolume, totalSurfaceArea, avSmoothness);
            Debug.Log("totalSurfaceArea - " + totalSurfaceArea);
            Debug.Log("totalVolume - " + totalVolume);
        }
    }
    private static void getCentrePoint(AudioReverbZone currentRZ)
    {
        Transform[] childrenTransforms = currentRZ.GetComponentsInChildren<Transform>();

        foreach (var item in childrenTransforms)
        {
            if (item.CompareTag("CentrePoint"))
            {
                m_CurrentCentrePoint = item.gameObject;
            }
        }
    }
    private static List<Vector3> getPoints(int samples)
    {
        List<Vector3> points = new List<Vector3>(samples);
        float phi = Mathf.PI * (3.0f - Mathf.Sqrt(5));  //calculates golden angle in radians
        float fSamples = System.Convert.ToSingle(samples);
        for (int i = 0; i < samples; i++)
        {
            float fI = System.Convert.ToSingle(i);
            float y = 1.0f - (fI / (fSamples - 1.0f)) * 2.0f;  //  y goes from 1 to -1
            float radius = Mathf.Sqrt(1.0f - y * y);   //  radius at y
            float theta = phi * i;  //golden angle increment
            float x = Mathf.Cos(theta) * radius;
            float z = Mathf.Sin(theta) * radius;

            Vector3 pnt = new Vector3(x, y, z);

            points.Add(pnt);
        }

        //Debug.Log(points.Count);
        return points;
    }
    private static IList<RaycastHit> raycasts(IList<Vector3> points)
    {
        Vector3 origin = m_CurrentCentrePoint.transform.position;
        IList<RaycastHit> collisions = new List<RaycastHit>(points.Count);

        int i = 0;
        foreach (var direction in points)
        {
            Ray ray = new Ray(origin, direction);
            //Debug.DrawRay(origin, direction, Color.red);
            RaycastHit collision = new RaycastHit();

            if (Physics.Raycast(ray, out collision, Mathf.Infinity))
            {
                //Debug.Log("hit" + collision.point);
                collisions.Add(collision);
            }
            else
            {
                //Debug.DrawLine(ray.origin, collision.point, Color.red);
                //Debug.DrawRay(origin, direction, Color.red);
                collisions.Add(new RaycastHit());
                //System.Diagnostics.Debugger.Break();
            }
            i++;
        }
        return collisions;
    }
    private static List<int> DelaunayTriangulation(List<Vector3> points)
    {
        IList<Vector2> points2D = pointsTo2D(points);
        List<int> triangles = new List<int>(points.Count);

        Delaunator triangulator = new Delaunator(DelaunatorSharp.Unity.Extensions.DelaunatorExtensions.ToPoints(points2D));
        triangles.AddRange(triangulator.Triangles);

        IPoint[] hullPoints = triangulator.GetHullPoints();
        int hullCount = hullPoints.Length;
        int firstIndex = 0;
        int prevIndex = 0;
        for (int i = 0; i < hullCount; i++)
        {
            Vector2 hullPoint = new Vector2(System.Convert.ToSingle(hullPoints[i].X), System.Convert.ToSingle(hullPoints[i].Y));
            int hullIndex = points2D.IndexOf(hullPoint);
            if (prevIndex != 0)
            {
                if (hullIndex != -1)
                {
                    triangles.Add(prevIndex);
                    triangles.Add(points.Count - 1);
                    triangles.Add(hullIndex);
                }
            }
            if (i == hullCount - 1)
            {
                triangles.Add(hullIndex);
                triangles.Add(points.Count - 1);
                triangles.Add(firstIndex);
            }

            if (firstIndex == 0)
                firstIndex = hullIndex;
            prevIndex = hullIndex;
        }

        return triangles;
    }
    private static List<Vector2> pointsTo2D(IList<Vector3> points)
    {
        List<Vector2> points2D = new List<Vector2>(points.Count);
        int badPoints = 0;
        foreach (var point in points)
        {
            if (point.x == 0.0f && point.z == 0.0f && point.y < 0.0f)
            {
                badPoints++;
                continue;
            }

            float x = point.x;
            float z = point.z;

            if (point.y < 0.0f)
            {
                float a = 1 / (Mathf.Pow(point.x, 2) + Mathf.Pow(point.z, 2));
                x = a * x;
                z = a * z;
            }

            points2D.Add(new Vector2(x, z));
        }
        return points2D;
    }
    private static List<Vector3> rcCollisionPoints(IList<RaycastHit> raycasts)
    {
        List<Vector3> collisionPoints = new List<Vector3>(raycasts.Count);

        foreach (var item in raycasts)
        {
            collisionPoints.Add(item.point);
        }

        return collisionPoints;
    }
    private static void geometricCalculations(List<int> triangles, List<Vector3> collisionPoints, IList<RaycastHit> raycastCollisions,out float[] avAbsorptionSpec, out float totalVolume, out float totalSurfaceArea, out float avSmoothness)
    {
        totalVolume = 0;
        totalSurfaceArea = 0;
        float[] absorptionSpectrum;
        avAbsorptionSpec = new float[7];
        float smoothness;
        avSmoothness = 0;

        for (int i = 0; i < triangles.Count; i = i + 3)
        {
            Vector3 point1 = collisionPoints[triangles[i]];
            Vector3 point2 = collisionPoints[triangles[i + 1]];
            Vector3 point3 = collisionPoints[triangles[i + 2]];

            float surfaceArea = calcSurfaceArea(point1, point2, point3);
            totalSurfaceArea += surfaceArea;

            triangleAbsorptionSpectrum(i, triangles, raycastCollisions, surfaceArea, out absorptionSpectrum, out smoothness);
            avAbsorptionSpec = addAbsorptionSpectrums(avAbsorptionSpec, absorptionSpectrum);
            avSmoothness += smoothness;

            totalVolume += calcVolume(point1, point2, point3, surfaceArea);
        }

        return;
    }
    private static float calcSurfaceArea(Vector3 point1, Vector3 point2, Vector3 point3)
    {
        //heron's formula / 0.1 is scale factopr for coordinates
        float side1 = SCALEFACTOR * Vector3.Distance(point1, point2); // distance round triangle
        float side2 = SCALEFACTOR * Vector3.Distance(point2, point3);
        float side3 = SCALEFACTOR * Vector3.Distance(point3, point1);

        float s = (side1 + side2 + side3) / 2;
        float surfaceArea = Mathf.Sqrt(s * ((s - side1) * (s - side2) * (s - side3)));

        //Debug.Log("Surface Area: " + surfaceArea);

        return surfaceArea;
    }
    private static float calcVolume(Vector3 point1, Vector3 point2, Vector3 point3, float surfaceArea)
    {
        float volume;
        Vector3 origin = m_CurrentCentrePoint.transform.position;

        volume = Mathf.Abs((Vector3.Dot(point1 - origin, Vector3.Cross(point2 - origin, point3 - origin))) / 6); //volume of a tetrahedrom
        volume = volume * Mathf.Pow(SCALEFACTOR, 3); //scaling

        return volume;
    }
    private static void triangleAbsorptionSpectrum(int i, IList<int> triangles, IList<RaycastHit> raycastCollisions, float surfaceArea, out float[] absorptionSpectrum, out float smoothness)
    {
        absorptionSpectrum = new float[7];
        smoothness = 0;

        // cycles through three indecies
        for (int point = i; point < (i + 3); point++) 
        {
            //gets collider ray interacts with
            Collider collider = raycastCollisions[triangles[i]].collider;
            if(collider != null)
            {
                // gets gameobject's material
                GameObject gameobjectAtPoint = collider.gameObject; 
                MeshRenderer objectMeshRenderer = gameobjectAtPoint.GetComponent<MeshRenderer>();
                Material objectMaterial = objectMeshRenderer.sharedMaterial;
                if (objectMaterial != null) //if material look up material name
                {
                    //sums absorption spectrums of all three points
                    float[] pointAbsorptionSpec = absorptionSpectrumLookup(objectMaterial); 
                    absorptionSpectrum = addAbsorptionSpectrums(absorptionSpectrum, pointAbsorptionSpec);
                    //adds smoothness of points
                    smoothness += objectMaterial.GetFloat("_Glossiness");
                }
            }
        }
        for (int index = 0; index < 7; index++)
        {
            //divides values by three
            absorptionSpectrum[index] = (absorptionSpectrum[index] / 3) * surfaceArea;
            //Debug.Log(absorptionSpectrum[index]);
        }

        smoothness = (smoothness / 3) * surfaceArea;
        //Debug.Log(smoothness);
        return;
    }
    static int counter = 0;
    private static float[] absorptionSpectrumLookup(Material material) //looks up material absorption spec
    {
        float[] absorptionSpectrum;

        //Debug.Log("spike " + (++counter).ToString() + " " + meterial.name);

        if (material.name.Contains("metal"))
        {
            absorptionSpectrum = new float[7] { 0.01f, 0.01f, 0.01f, 0.02f, 0.02f, 0.02f, 0.02f };
        }
        else if (material.name.Contains("wood"))
        {
//            Debug.Log("wood");
            absorptionSpectrum = new float[7] { 0.27f, 0.23f, 0.22f, 0.15f, 0.10f, 0.07f, 0.06f };
        }
        else if (material.name.Contains("brick"))
        {
//            Debug.Log("brick");
            absorptionSpectrum = new float[7] { 0.03f, 0.03f, 0.03f, 0.04f, 0.05f, 0.07f, 0.07f };
        }
        else if (material.name.Contains("plaster"))
        {
//            Debug.Log("plaster");
            absorptionSpectrum = new float[7] { 0.15f, 0.10f, 0.06f, 0.04f, 0.04f, 0.05f, 0.05f };
        }
        else if (material.name.Contains("tile"))
        {
//            Debug.Log("tile");
            absorptionSpectrum = new float[7] { 0.01f, 0.01f, 0.02f, 0.02f, 0.02f, 0.02f, 0.02f };
        }
        else if (material.name.Contains("absorption"))
        {
//            Debug.Log("absorption");
            absorptionSpectrum = new float[7] { 0.22f, 0.60f, 0.92f, 0.90f, 0.88f, 0.88f, 0.88f };
        }
        else if (material.name.Contains("concrete"))
        {
//            Debug.Log("concrete");
            absorptionSpectrum = new float[7] { 0.01f, 0.01f, 0.02f, 0.02f, 0.02f, 0.05f, 0.05f };
        }
        else if (material.name.Contains("glass"))
        {
//            Debug.Log("glass");
            absorptionSpectrum = new float[7] { 0.08f, 0.04f, 0.03f, 0.03f, 0.02f, 0.02f, 0.02f };
        }
        else if (material.name.Contains("carpet"))
        {
            //            Debug.Log("carpet");
            absorptionSpectrum = new float[7] { 0.07f, 0.31f, 0.49f, 0.81f, 0.66f, 0.54f, 0.48f };
        }
        else if (material.name.Contains("curtain"))
        {
            //            Debug.Log("curtain");
            absorptionSpectrum = new float[7] { 0.3f, 0.45f, 0.65f, 0.56f, 0.59f, 0.71f, 0.71f };
        }
        else
        {
//            Debug.Log("unknown");
            absorptionSpectrum = new float[7] { 0.3f, 0.3f, 0.3f, 0.3f, 0.3f, 0.3f, 0.3f };
        }
        return absorptionSpectrum;
    }
    private static void setParameters(AudioReverbZone reverbZone, float[] avAbsorptionSpec, float totalVolume, float totalSurfaceArea, float avSmoothness)
    {
        //average absorption at frequency bands including air absorption.
        double highFreq = (((avAbsorptionSpec[0] / totalSurfaceArea) + (avAbsorptionSpec[1] / totalSurfaceArea)) / 2) + (0.002 * totalVolume); 
        double midFreq = (((avAbsorptionSpec[2] / totalSurfaceArea) + (avAbsorptionSpec[3] / totalSurfaceArea) + (avAbsorptionSpec[4] / totalSurfaceArea)) / 3);
        double lowFreq = (((avAbsorptionSpec[5] / totalSurfaceArea) + (avAbsorptionSpec[6] / totalSurfaceArea)) / 2);

        //Room - room effect level at mid frequencies (-10000 - 0)
        reverbZone.room = System.Convert.ToInt32(0 - (midFreq * 10000));
        reverbZone.roomHF = System.Convert.ToInt32(0 - (highFreq * 10000));
        reverbZone.roomLF = System.Convert.ToInt32(0 - (lowFreq * 10000));
        
        //decay time - sabines equation
        double decayTime = (0.161 * totalVolume) / (avAbsorptionSpec[3] + (4 * 0.003 * totalVolume));
        reverbZone.decayTime = System.Convert.ToSingle(decayTime);

        //Decay HF Ratio
        double decayTimeMid = (0.161 * totalVolume) / (midFreq * totalSurfaceArea);
        double decayTimeHigh = (0.161 * totalVolume) / (highFreq * totalSurfaceArea);
        reverbZone.decayHFRatio = System.Convert.ToSingle(decayTimeHigh / decayTimeMid);

        //Reflections - scale 1000 - -10000
        reverbZone.reflections = System.Convert.ToInt32(0 - (midFreq * 11000));

        //reflections delay
        reverbZone.reflectionsDelay = System.Convert.ToSingle( 2 * ((reverbZone.minDistance * SCALEFACTOR) / 343)); //delay = distance / speed of sound

        //Reverb
        reverbZone.reverb = System.Convert.ToInt32(1000 - (midFreq * 12000) );

        //ReverbDelay
        reverbZone.reverbDelay = System.Convert.ToSingle((decayTime / 2) / 100); // check

        //HF LF Reference
        reverbZone.LFReference = 250;
        reverbZone.HFReference = 4000;

        //Diffusion
        reverbZone.diffusion = 100 - ((avSmoothness / totalSurfaceArea) * 30) ; //mapped between 70 - 100

        //Density
        reverbZone.density = 100;

        //preset
        reverbZone.reverbPreset = AudioReverbPreset.User;
    }
    private static float[] addAbsorptionSpectrums(float[] x, float[] y)
    {
        float[] output = new float[7];

        for (int i = 0; i < 7; i++)
        {
            //Debug.Log(x[i]);
            //Debug.Log(y[i]);
            output[i] = x[i] + y[i];
            //Debug.Log(y[i]);
        }
        return output;
    }
}


