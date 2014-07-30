package org.liufeng.course.service;

import java.net.*;
import java.io.*;

public class MessageRouter {
	public static final String HOST = "127.0.0.1";
	public static final int PORT = 9999;
	
	public static void sendMessage(String msg) throws UnknownHostException, IOException {
		Socket socket = new Socket(HOST, PORT);
		PrintWriter os = new PrintWriter(socket.getOutputStream());
		os.println(msg);
		socket.close();
	}
}
