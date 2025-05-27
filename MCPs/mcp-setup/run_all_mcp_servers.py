import asyncio
import os
import sys
import glob

# Path to the directory containing this script and the server scripts
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

async def log_output_continuously(stream, prefix):
    """Continuously logs output from a stream until EOF or cancellation."""
    try:
        while True:
            line = await stream.readline()
            if line:
                print(f"     [{prefix}] {line.decode(errors='ignore').strip()}")
            else:
                # End of stream
                print(f"     [{prefix}] Stream closed (EOF).")
                break
    except asyncio.CancelledError:
        print(f"     [{prefix}] Logging task cancelled.")
        # Do not re-raise, allow the task to be considered "handled" upon cancellation.
    except Exception as e:
        print(f"     [{prefix}] Error reading output: {type(e).__name__} - {e}")
    # finally:
        # print(f"     [{prefix}] Logging task finished execution.") # Optional: for debugging

async def start_server(server_script_path):
    server_name = os.path.basename(server_script_path)
    # This print is moved to main after successful launch confirmation
    # print(f"Attempting to start {server_name}...") 
    
    env = os.environ.copy()
    env["MCP_TRANSPORT"] = "sse"
    
    python_path = env.get("PYTHONPATH", "")
    script_dir_abs_path = os.path.abspath(SCRIPT_DIR)
    if script_dir_abs_path not in [os.path.abspath(p) for p in python_path.split(os.pathsep) if p]:
        env["PYTHONPATH"] = f"{script_dir_abs_path}{os.pathsep}{python_path}".strip(os.pathsep)

    process = await asyncio.create_subprocess_exec(
        sys.executable, server_script_path,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=SCRIPT_DIR,
        env=env
    )
    
    # Confirmation of launch will be handled in main based on result of this coroutine
    # print(f"  -> Launched {server_name} (PID: {process.pid}). It will run with MCP_TRANSPORT=sse.")
    
    stdout_log_task = asyncio.create_task(log_output_continuously(process.stdout, f"{server_name} STDOUT"))
    stderr_log_task = asyncio.create_task(log_output_continuously(process.stderr, f"{server_name} STDERR"))
    
    return process, stdout_log_task, stderr_log_task, server_name

async def main():
    print("Starting all MCP servers asynchronously with MCP_TRANSPORT=sse...\\n")
    
    server_pattern = os.path.join(SCRIPT_DIR, "*_mcp_server.py")
    server_scripts = glob.glob(server_pattern)
    
    current_script_name = os.path.basename(__file__)
    server_scripts = [s for s in server_scripts if os.path.basename(s) != current_script_name]

    if not server_scripts:
        print(f"No MCP server scripts found matching '{server_pattern}' in {SCRIPT_DIR}.")
        print("Ensure server scripts end with '_mcp_server.py'.")
        return
        
    print(f"Found the following server scripts to start: {', '.join([os.path.basename(s) for s in server_scripts])}\\n")
    
    start_coroutines = [start_server(script_path) for script_path in server_scripts]
    
    results = await asyncio.gather(*start_coroutines, return_exceptions=True)
    
    launched_servers_info = [] # Stores (process, stdout_task, stderr_task, server_name)

    for i, res in enumerate(results):
        # Fallback server name if start_server failed before returning server_name
        script_path_for_fallback = server_scripts[i]
        server_name_fallback = os.path.basename(script_path_for_fallback)

        if isinstance(res, tuple) and len(res) == 4 and isinstance(res[0], asyncio.subprocess.Process):
            process, stdout_task, stderr_task, server_name = res
            launched_servers_info.append((process, stdout_task, stderr_task, server_name))
            print(f"  -> Successfully started and monitoring: {server_name} (PID: {process.pid})")
        elif isinstance(res, Exception):
            print(f"  !!! Failed to start {server_name_fallback}: {res} !!!")
        else: # Should not happen if start_server is structured correctly
            print(f"  !!! Unknown issue starting {server_name_fallback}: {res} !!!")

    if not launched_servers_info:
        print("\\nNo servers were successfully started. Exiting.")
        return

    all_logging_tasks = []
    for _, stdout_task, stderr_task, _ in launched_servers_info:
        all_logging_tasks.append(stdout_task)
        all_logging_tasks.append(stderr_task)

    print("\\n--------------------------------------------------------------------")
    print(f"{len(launched_servers_info)} server(s) are running. Output is being logged below.")
    print("Press Ctrl+C to stop all servers and exit.")
    print("--------------------------------------------------------------------\\n")

    try:
        # Keep main alive while logging tasks run.
        # This will complete when all logging tasks complete (e.g., all servers exit and close pipes)
        # or if this gather call itself is cancelled (e.g., due to KeyboardInterrupt on main coroutine).
        if all_logging_tasks:
            await asyncio.gather(*all_logging_tasks, return_exceptions=True) # return_exceptions to log errors from loggers
        print("\\nAll server logging tasks have completed (servers may have exited).")

    except KeyboardInterrupt:
        print("\\n\\nCtrl+C received. Initiating shutdown of all servers...")
        for process, _, _, server_name in launched_servers_info:
            if process.returncode is None: # If process is still running
                print(f"  Terminating {server_name} (PID: {process.pid})...")
                try:
                    process.terminate()
                except Exception as e_term:
                    print(f"    Error during terminate for {server_name}: {e_term}")
        
        print("  Waiting for server processes to exit...")
        wait_for_exit_tasks = []
        for process, _, _, server_name in launched_servers_info:
            if process.returncode is None:
                async def wait_and_log_exit(p, s_name):
                    try:
                        ret_code = await p.wait()
                        print(f"    Server {s_name} (PID: {p.pid}) exited with code {ret_code}.")
                    except Exception as e_wait:
                        print(f"    Error waiting for {s_name} (PID: {p.pid}) to exit: {e_wait}")
                wait_for_exit_tasks.append(wait_and_log_exit(process, server_name))
        
        if wait_for_exit_tasks:
            await asyncio.gather(*wait_for_exit_tasks, return_exceptions=True) # Errors handled in wait_and_log_exit

        print("  All server processes signaled for termination and awaited.")

    except Exception as e_main:
        print(f"\\nUnexpected error in main monitoring loop: {e_main}")
        # Potentially trigger shutdown here as well if needed, though finally should cover cleanup.

    finally:
        print("\\n--------------------------------------------------------------------")
        print("Performing final cleanup...")

        # Cancel any logging tasks that might still be technically active
        # (e.g., if a stream read was stuck and didn't respond to earlier EOF/cancellation)
        active_logging_tasks = [task for task in all_logging_tasks if not task.done()]
        if active_logging_tasks:
            print(f"  Cancelling {len(active_logging_tasks)} outstanding logging task(s)...")
            for task in active_logging_tasks:
                task.cancel()
            # Wait for these cancellations to be processed
            await asyncio.gather(*active_logging_tasks, return_exceptions=True)
            print("  Outstanding logging tasks cancellation processed.")
        #else:
            #print("  All logging tasks were already completed or handled.")

        # Clean up subprocess transports
        if launched_servers_info:
            print("  Cleaning up server process transports...")
            for proc, _, _, server_name in launched_servers_info:
                proc_identifier = f"{server_name} (PID: {proc.pid if proc else 'Unknown'})"
                
                # Ensure process is truly stopped if it wasn't by terminate/wait
                if proc and proc.returncode is None:
                    print(f"     Warning: Process {proc_identifier} still running during final cleanup. Attempting to kill.")
                    try:
                        proc.kill() # Force kill
                        await asyncio.wait_for(proc.wait(), timeout=5.0) # Give it a moment to die
                        print(f"     Process {proc_identifier} killed, exit code {proc.returncode}.")
                    except asyncio.TimeoutError:
                        print(f"     Warning: Timeout waiting for {proc_identifier} to die after kill.")
                    except Exception as e_kill:
                        print(f"     Warning: Error killing process {proc_identifier}: {e_kill}")

                # Close stdin if it exists and is open
                if proc and proc.stdin and hasattr(proc.stdin, 'close') and not proc.stdin.is_closing():
                    try:
                        proc.stdin.close()
                    except Exception as e:
                        print(f"     Warning: Error closing stdin for {proc_identifier}: {e}")
                
                # Close the transport
                if proc and hasattr(proc, '_transport') and proc._transport:
                    if hasattr(proc._transport, 'is_closing') and not proc._transport.is_closing():
                        try:
                            proc._transport.close()
                        except Exception as e:
                            print(f"     Warning: Error closing transport for {proc_identifier}: {e}")
                    elif not hasattr(proc._transport, 'is_closing'): # Fallback for other transport types/versions
                         try:
                            proc._transport.close()
                         except Exception as e: # pylint: disable=broad-except
                            print(f"     Warning: Error closing transport (no is_closing check) for {proc_identifier}: {e}")
                # else:
                    # print(f"     Info: No transport to close or already closed for {proc_identifier}.")


            print("  Allowing a brief moment for transport cleanup tasks to finalize...")
            await asyncio.sleep(0.25) # Slightly longer for more complex cleanup

        print("Cleanup complete. Runner script will now exit.")
        print("--------------------------------------------------------------------")

if __name__ == "__main__":
    if sys.platform == "win32" and isinstance(asyncio.get_event_loop_policy(), asyncio.WindowsSelectorEventLoopPolicy):
        asyncio.set_event_loop_policy(asyncio.WindowsProactorEventLoopPolicy())
    
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        # This secondary KeyboardInterrupt catch is for asyncio.run(main()) itself.
        # The main() coroutine has its own internal KeyboardInterrupt handling.
        # This ensures that if Ctrl+C happens *during* the setup of asyncio.run or very early,
        # the script still acknowledges it.
        print("\\nRunner script execution directly interrupted. Exiting.")
    except Exception as e_global:
        print(f"Global unhandled exception: {e_global}")
        # Consider logging this to a file or more robustly for production scripts
