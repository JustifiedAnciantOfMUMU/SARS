using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;
using DelaunatorSharp;

public class CreateRoomMeshV2
{
    private static AudioReverbZone[] m_ReverbZones;
    private static GameObject m_CurrentCentrePoint;
    private const float SCALEFACTOR = 0.1f;

    [MenuItem("Tools/Geometry/CreateRoomMeshV2")]
    private static void main()
    {
        m_ReverbZones = Object.FindObjectsOfType<AudioReverbZone>();

        foreach (var item in m_ReverbZones)
        {
            getCentrePoint(item);
            List<Vector3> points = getPoints(65535); // num points|| max num = 65535 due to int variable type
            // IList<Vector2> points2D = pointsTo2D(points);
            List<Vector3> drape;
            int[] allTriangles = calculatedrapeMesh(points, out drape); //DelaunayTriangulation(points);

            IList<RaycastHit> raycastCollisions = raycasts(points);

            Mesh mesh = new Mesh();
            List<Vector3> CollisionPoints = rcCollisionPoints(raycastCollisions);
            mesh.SetVertices(CollisionPoints);
            mesh.SetTriangles(allTriangles, 0);
            m_CurrentCentrePoint.GetComponent<MeshFilter>().mesh = mesh;

            List<int> listAllTriangles = new List<int>(allTriangles);
            float SA = geometricCalculations(listAllTriangles, CollisionPoints, raycastCollisions);
        }
    }

    private static void getCentrePoint(AudioReverbZone currentRZ)
    {
        Transform[] childrenTransforms = currentRZ.GetComponentsInChildren<Transform>();

        foreach (var item in childrenTransforms)
        {
            if (item.CompareTag("CentrePoint"))
            {
                //removing componenets
                if (item.gameObject.GetComponent<MeshRenderer>())
                {
                    //GameObject.DestroyImmediate(item.gameObject.GetComponent<MeshRenderer>());
                }
                if (item.gameObject.GetComponent<MeshFilter>())
                {
                    //GameObject.DestroyImmediate(item.gameObject.GetComponent<MeshFilter>());
                }
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

            //if (y < 0 && x == z && x == 0)
            //    System.Diagnostics.Debugger.Break();

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
            //Debug.Log(direction);
            Ray ray = new Ray(origin, direction);
            Debug.DrawRay(origin, direction, Color.red);
            RaycastHit collision = new RaycastHit();

            if (Physics.Raycast(ray, out collision, Mathf.Infinity))
            {
                Debug.Log("hit" + collision.point);
                Debug.DrawLine(ray.origin, collision.point, Color.red);
                collisions.Add(collision);
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

    private static int[] calculatedrapeMesh(List<Vector3> points, out List<Vector3> meshPoints)
    {
        meshPoints = new List<Vector3>(points.Count);
        IList<Vector2> points2D = pointsTo2D(points);

        Delaunator triangulator = new Delaunator(DelaunatorSharp.Unity.Extensions.DelaunatorExtensions.ToPoints(points2D));

        for (int i = 0; i < points2D.Count; i++)
        {
            meshPoints.Add(new Vector3(points2D[i].x, points[i].y, points2D[i].y));
        }


        return triangulator.Triangles;
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

    private static float geometricCalculations(List<int> triangles, List<Vector3> collisionPoints, IList<RaycastHit> raycastCollisions)
    {
        float totalSurfaceArea = 0.0f;
        float totalVolume = 0.0f;
        float[] AvAbsorptionSpec = new float[7];
        float alpha = 0.0f;

        for (int i = 0; i < triangles.Count; i = i + 3)
        {
            Vector3 point1 = collisionPoints[triangles[i]];
            Vector3 point2 = collisionPoints[triangles[i + 1]];
            Vector3 point3 = collisionPoints[triangles[i + 2]];

            float surfaceArea = calcSurfaceArea(point1, point2, point3);
            //float[] absorptionSpectrum = objectAbsorptionSpectrum(i, triangles, raycastCollisions, surfaceArea);
            //AvAbsorptionSpec  = absorptionSpectrum;
            //alpha = absorptionSpectrum[3] * surfaceArea;
            totalVolume += calcVolume(point1, point2, point3, surfaceArea);
            totalSurfaceArea += surfaceArea;
        }
        Debug.Log("SF: " + totalSurfaceArea);
        Debug.Log("TV: " + totalVolume);

        return totalSurfaceArea;
    }

    private static float calcSurfaceArea(Vector3 point1, Vector3 point2, Vector3 point3)
    {
        //heron's formula / 0.1 is scale factopr for coordinates
        float side1 = SCALEFACTOR * Vector3.Distance(point1, point2); // distance round triangle
        float side2 = SCALEFACTOR * Vector3.Distance(point2, point3);
        float side3 = SCALEFACTOR * Vector3.Distance(point3, point1);

        float s = (side1 + side2 + side3) / 2;
        float surfaceArea = Mathf.Sqrt(s * ((s - side1) * (s - side2) * (s - side3)));

        Debug.Log("Surface Area: " + surfaceArea);

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
}

