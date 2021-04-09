using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;
using DelaunatorSharp;

public class DrapeVislisation
{
    private static AudioReverbZone[] m_ReverbZones;
    private static GameObject m_CurrentCentrePoint;

    [MenuItem("Tools/Geometry/Drape Visulisation")]
    private static void main()
    {
        m_ReverbZones = Object.FindObjectsOfType<AudioReverbZone>();

        foreach (var item in m_ReverbZones)
        {
            getCentrePoint(item);
            List<Vector3> points = getPoints(250); // num points|| max num = 65535 due to int variable type
            List<Vector3> drape;
            int[] allTriangles = calculatedrapeMesh(points, out drape); //DelaunayTriangulation(points);

            IList<RaycastHit> raycastCollisions = raycasts(points);

            Mesh mesh = new Mesh();
            mesh.SetVertices(drape);
            mesh.SetTriangles(allTriangles, 0);
            m_CurrentCentrePoint.GetComponent<MeshFilter>().mesh = mesh;
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

        Debug.Log(points.Count);
        return points;
    }
    private static IList<RaycastHit> raycasts(IList<Vector3> points)
    {
        Vector3 origin = m_CurrentCentrePoint.transform.position;
        IList<RaycastHit> collisions = new List<RaycastHit>(points.Count);

        int i = 0;
        foreach (var direction in points)
        {
            Debug.Log(direction);
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
}


