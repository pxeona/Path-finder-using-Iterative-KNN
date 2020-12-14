package rmit;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.sql.Connection;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import db.Dataset;
import spatialindex.spatialindex.*;
import spatialindex.rtree.*;

public class IKNN {
	private RTree tree;
	private TreeMap<Double,Point> T;
	private HashMap<String,String> trips;
	private Connection conn = null; 
	private Dataset ds = null;
	public long iotime;
	public long candis;
	
	// compare the distance of two points to given point
	private double compareDistance(Point a, Point b, Point c){
		double dist_ab = a.getMinimumDistance(b);
		double dist_ac = a.getMinimumDistance(c);
		
		if(dist_ab > dist_ac){
			return 1;
		}
		return 0;		
	}
	
	public IKNN(RTree tree, String t, Dataset d, Connection c) throws SQLException{
		// Initializing values
		this.tree = tree;
		this.trips = this.readTrips(t);
		this.ds = d;
		this.conn = c;
	}
	private HashMap<String, String> readTrips(String FilePath){
	  // Hash to store POI and corresponding trajectory IDs
	  HashMap<String, String> trips = new HashMap<String, String>();
	  try {
	   InputStreamReader read = new InputStreamReader(new FileInputStream(FilePath), "utf-8");
	   BufferedReader reader = new BufferedReader(read);
	   String line;
	   String []arr = null;
	   String tmp = null;
	   while ((line = reader.readLine()) != null) {
		arr = line.split(",");
	    // Parsing value and adding into hash map
	    tmp = Double.parseDouble(arr[1])+","+Double.parseDouble(arr[2]);

	    if(trips.containsKey(tmp)){
	    	if(!trips.get(tmp).contains(arr[0])){
	    		trips.put(tmp, trips.get(tmp)+","+arr[0]);
	    	}    	
	    }
	    else{
	    	trips.put(tmp, arr[0]);
	    }
	    
	   } reader.close();
	  } catch (Exception e) {
	   e.printStackTrace();
	  }
	  return trips;
	}
		
	public String computeIKNN(Point []points, int k) throws SQLException{
		
		iotime = 0;
		String output = "";
		int lambda = k;
		
		// stores candidate trajectories
		HashMap<String, ArrayList<Point>> candidates = new HashMap<String, ArrayList<Point>>();
		// Stores upper bound values
		ArrayList<Point> UB_points = new ArrayList<Point>();
		for (int i = 0; i < points.length; i++) {
			UB_points.add(null);
		}
		int index = 0;
		int iteration=1;
		int counter = 0;
		int check = 0;
		
		System.out.println("Number of Query Points:" + points.length);
		

		//filter the impossible trajectories.
		while (check == 0){
			System.out.println("Iteration: " + iteration);
			for (int i = 0; i < points.length; i++) {
				// Starts search from first query point and iterates
				index = i;
				ArrayList<Point> S = new ArrayList<Point>();
				long startTime = System.currentTimeMillis();
				// finds lambda nearest points to given query point
				this.getIntersectingPoints(S, lambda, points[i]);
				//System.out.println("Result Point: " + S.size());
				long stopTime = System.currentTimeMillis();
				long elapsedTime = stopTime - startTime;
				iotime += elapsedTime;
				counter = S.size();
				System.out.println("Query Point Index: " + index + " Points: " + S.size());
				String trajectPassingthruQueryLoc = null;
				
				for (Point point : S) 
				{	
					trajectPassingthruQueryLoc = null;
					if(UB_points.get(i) == null)
					{
						UB_points.set(i, point);
					}
					else if(points[i].getMinimumDistance(UB_points.get(i)) < points[i].getMinimumDistance(point))
					{
						UB_points.set(i, point);
					}
					
					String queryPoint[] = point.toString().split(":");
					trajectPassingthruQueryLoc = trips.get(Double.parseDouble(queryPoint[0]) + "," + Double.parseDouble(queryPoint[1]));
					if(trajectPassingthruQueryLoc != null)
					{
						String trajIds[] = trajectPassingthruQueryLoc.split(",");
						for(int j = 0; j < trajIds.length ; j++)
						{
						    if(candidates.containsKey(trajIds[j]))
						    {
						    	if(!candidates.get(trajIds[j]).contains(point))
						    	{
						    		candidates.get(trajIds[j]).add(point);
						    	}
						    }
						    else
						    {
						    	ArrayList<Point> tracedPoint = new ArrayList<Point>();
						    	tracedPoint.add(point);
						    	candidates.put(trajIds[j],tracedPoint);
						    }
						}
					}

					//2. finds trajectory ids for each point and stores trajectories in candidate set, i.e., update candidates.
				
				}
			}
			
			if(candidates.size() >= k){
				
				PriorityQueue<Double> LB = new PriorityQueue<>(Collections.reverseOrder());
				// add candidate trajectory distance to query in descending order
				for (Map.Entry<String, ArrayList<Point>> entry : candidates.entrySet()) {
					double tmp_dist = computeLowerBound(entry.getValue(), points);
					LB.add(tmp_dist);
				}
				
				for (int i = 1; i < k; i++) {
					LB.poll();
				}
				// Choose k-th lower bound value
				double k_LB = LB.peek();
				
				// computes upper bound value using farthest point we found so far for each query point.
				double UB = computeUpperBound(UB_points, points);
				
				System.out.println("Current UB: " + UB + " Current k-th LB: " + k_LB );
				
				if(k_LB >= UB){
					check = 1;
				}
				
			}
			// increase lambda value by 50 and continues
			//System.out.println("Candidates:" + candidates.size());
			lambda += 50;
			iteration++;
			System.out.println("Candidates:" + candidates.size() + "\n");
		}
		
		System.out.println("End of candidate generation");
		

		//refine the candidate trajectories
		
		PriorityQueue<Candidate> resultSet = new PriorityQueue<>();// stores top-k results		
		PriorityQueue<Candidate> sorted_candidates = new PriorityQueue<>(Collections.reverseOrder());// candidates are stored in descending order of upper bound distance
		for (Map.Entry<String, ArrayList<Point>> entry : candidates.entrySet()) {
			double candidate_ub = computeCandidateUpperBound(entry.getValue(), points, UB_points);
			sorted_candidates.add(new Candidate(candidate_ub, entry.getKey()));
		}
		int c = 0;
		while( sorted_candidates.peek() != null) 
		{
			String trajId = sorted_candidates.peek().getID();
		    double kBCTMin;
		    double sim = computeCandidateDistance(trajId, points);
		    if(resultSet.size() < k)
		    {
		    	resultSet.add(new Candidate(sim, trajId));
		    	sorted_candidates.poll();
		    }
		    else
		    {
		    	kBCTMin = resultSet.peek().getDistance();
		    	if(sim > kBCTMin)
		    	{
		    		resultSet.poll();
		    		resultSet.add(new Candidate(sim, trajId));
		    	}
		    	
	    		sorted_candidates.poll();
	    		
	    		if(sorted_candidates.peek() != null)
	    		{
			    	if(kBCTMin >= sorted_candidates.peek().getDistance())
			    	{
			    		break;
			    	}
	    		}
		    }
		    
		}		

		// format top-k results for printing
		while( resultSet.peek() != null) {
			Candidate can = resultSet.poll();
			output += can.getID()+"\t"+can.getDistance()+"\n";
		}
		this.candis = candidates.size();
		//System.out.println("Iteration: "+iteration+ " Candidates: "+candidates.size() + " Points:" +counter);
		return output;
	}
	 
	public int getIntersectingPoints(ArrayList<Point> s, int lambda, Point p){
		 MyVisitor v =  new MyVisitor();
		 // finds the nearest lambda points for given query point p
		 tree.nearestNeighborQuery(lambda, p, v);
		 
		 for (Map.Entry<Integer, IShape> entry : v.answers.entrySet()) {
		    IShape value = entry.getValue();
		    double[] coord = value.getCenter();
		    s.add(new Point(coord));
		 }
		 return 1;
	} 
	
	// 4 compute the lower bound distance of k-th results
	private double computeLowerBound(ArrayList<Point> p, Point[] points)
	{
		double dist=0;
		double minDistBtwQueryAndTrajPoint;
		for(int i = 0; i < points.length; i++)
		{
			minDistBtwQueryAndTrajPoint = points[i].getMinimumDistance(p.get(0));
			for(Point point : p)
			{
				if(points[i].getMinimumDistance(point) > minDistBtwQueryAndTrajPoint)
				{
					continue;
				}
				else
				{
					minDistBtwQueryAndTrajPoint = points[i].getMinimumDistance(point);
				}
			}
			dist += 1/(Math.exp(minDistBtwQueryAndTrajPoint));
		}
		return dist;		
	}
	// 5 computes upper bound of unseen trajectories using corresponding points 
	private double computeUpperBound(ArrayList<Point> p, Point[] points)
	{
		double dist=0;
		double distBetweenQueryAndNN;
		for(int i = 0; i < points.length; i++)
		{
			distBetweenQueryAndNN = points[i].getMinimumDistance(p.get(i));
			dist += 1/(Math.exp(distBetweenQueryAndNN));
		}
		return dist;
	}
	// 6 compute the upper bound of candidates 
	private double computeCandidateUpperBound(ArrayList<Point> p, Point[] points, ArrayList<Point> UB_points)
	{
		double dist=0;
		double minDistBtwQueryAndTrajPoint;
		double distBetweenQueryAndNN;
		double distLB = 0;
		double distUB = 0;
		for(int i = 0; i < points.length; i++)
		{
			minDistBtwQueryAndTrajPoint = points[i].getMinimumDistance(p.get(0));
			for(Point point : p)
			{
				if(points[i].getMinimumDistance(point) > minDistBtwQueryAndTrajPoint)
				{
					continue;
				}
				else
				{
					minDistBtwQueryAndTrajPoint = points[i].getMinimumDistance(point);
				}
			}
			distLB += 1/(Math.exp(minDistBtwQueryAndTrajPoint));
			distBetweenQueryAndNN = points[i].getMinimumDistance(UB_points.get(i));
			distUB += 1/(Math.exp(distBetweenQueryAndNN));
		}
		dist = distLB + distUB;
		return dist;		
	}
	
	// 7 Compute actual distance of the trajectory to query points
	public double computeCandidateDistance(String id, Point[] points) throws SQLException
	{
		double dist=0;
		String arrKey[];
		String arrValue[];
		ArrayList<String> tripPoints = new ArrayList<String>();
		ArrayList<Point> tracedPoints = new ArrayList<Point>();
		double minDistBtwQueryAndTrajPoint;
		for (Map.Entry<String, String> entry : trips.entrySet())
		{
			tripPoints.clear();
			arrValue = entry.getValue().split(",");
			for(int i = 0; i < arrValue.length; i++)
			{
				tripPoints.add(arrValue[i]);
			}
			if(tripPoints.contains(id))
			{
				
				arrKey = entry.getKey().split(",");
				double[] coords = {Double.parseDouble(arrKey[0]),Double.parseDouble(arrKey[1])};
				Point point = new Point(coords);
				tracedPoints.add(point);
			}
		}
		
		for(int j = 0; j < points.length; j++)
		{

			minDistBtwQueryAndTrajPoint = points[j].getMinimumDistance(tracedPoints.get(0));
			for(Point tp : tracedPoints)
			{
				if(points[j].getMinimumDistance(tp) > minDistBtwQueryAndTrajPoint)
				{
					continue;
				}
				else
				{
					minDistBtwQueryAndTrajPoint = points[j].getMinimumDistance(tp);
				}
			}
			dist += 1/(Math.exp(minDistBtwQueryAndTrajPoint));
		}
		
		return dist;		
	}
}


