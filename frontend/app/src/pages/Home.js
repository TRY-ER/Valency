import React from 'react';
import { Link } from 'react-router-dom';
import { motion } from 'framer-motion';
import { useTheme } from '../contexts/ThemeContext';
import { fadeInUpVariants, fadeInRightVariants, fadeInLeftVariants } from '../components/animations/framerAnim';
import DnaViewer from '../components/DnaViewer';
import '../styles/Home.css';
import CanvasSample from '../components/sample_3d_viewer';

function Home() {
  const { theme } = useTheme();

  return (
    <div className="home-container">
      {/* Navigation Bar */}
      {/* Hero Section */}
      <section className="hero-section">
        <motion.div 
          className="hero-content"
          initial="hidden"
          animate="visible"
          variants={fadeInLeftVariants}
        >
          <h1 className="hero-title">Accelerate Drug Discovery with Your Intelligent AI Assistant</h1>
          <p className="hero-subtitle">
            Seamlessly automate protein analysis, identify drug candidates, and optimize targets. 
            Take control with intuitive chat-based tweaking and integrated manual tools.
          </p>
          <div className="hero-cta">
            <Link to="/demo" className="primary-btn">Request a Demo</Link>
            <Link to="/explore" className="secondary-btn">Explore the Assistant</Link>
          </div>
        </motion.div>
        
        <motion.div 
          className="hero-dna-model"
          initial="hidden"
          animate="visible"
          variants={fadeInRightVariants}
        >
          {/* <DnaViewer /> */}
          {/* <CanvasSample /> */}
        </motion.div>
      </section>

      {/* Why Valency Section */}
      <section className="why-section" id="why">
        <motion.h2 
          className="section-title"
          initial="hieden"
          whileInView="visible"
          viewport={{ once: true }}
          variants={fadeInUpVariants}
        >
          Why Valency?
        </motion.h2>
        <div className="problem-cards">
          <motion.div 
            className="problem-card glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={0}
          >
            <h3>Time-Consuming Processes</h3>
            <p>Traditional drug discovery involves labor-intensive manual processes that slow research progress.</p>
          </motion.div>
          
          <motion.div 
            className="problem-card glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={1}
          >
            <h3>Fragmented Tools</h3>
            <p>Researchers waste valuable time switching between disconnected tools and reconciling inconsistent data.</p>
          </motion.div>
          
          <motion.div 
            className="problem-card glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={2}
          >
            <h3>Steep Learning Curves</h3>
            <p>Complex software requires extensive training, creating barriers to effective utilization.</p>
          </motion.div>
        </div>
        <motion.div 
          className="solution-statement"
          initial="hidden"
          whileInView="visible"
          viewport={{ once: true }}
          variants={fadeInUpVariants}
        >
          <h2>We bridge the gap with an AI assistant that automates workflows while keeping you in control.</h2>
        </motion.div>
      </section>

      {/* How It Works Section */}
      <section className="how-works-section" id="how-it-works">
        <motion.h2 
          className="section-title"
          initial="hidden"
          whileInView="visible"
          viewport={{ once: true }}
          variants={fadeInUpVariants}
        >
          Meet Your AI Drug Discovery Co-Pilot
        </motion.h2>
        
        <div className="workflow-steps">
          <motion.div 
            className="workflow-step"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={0}
          >
            <div className="step-number">1</div>
            <div className="step-content">
              <h3>Connect & Configure</h3>
              <p>Seamlessly connect your data sources and configure your research parameters with natural language.</p>
            </div>
          </motion.div>
          
          <motion.div 
            className="workflow-step"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={1}
          >
            <div className="step-number">2</div>
            <div className="step-content">
              <h3>Automate & Analyze</h3>
              <p>Let the AI execute complex workflows across protein analysis, drug candidate identification, and target optimization.</p>
            </div>
          </motion.div>
          
          <motion.div 
            className="workflow-step"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={2}
          >
            <div className="step-number">3</div>
            <div className="step-content">
              <h3>Interact & Tweak via Chat</h3>
              <p>Guide the agent, adjust parameters, and request specific analyses through an intuitive chat interface.</p>
            </div>
          </motion.div>
          
          <motion.div 
            className="workflow-step"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={3}
          >
            <div className="step-number">4</div>
            <div className="step-content">
              <h3>Utilize Integrated Manual Tools</h3>
              <p>Access specialized tools directly within the platform for hands-on control when needed.</p>
            </div>
          </motion.div>
          
          <motion.div 
            className="workflow-step"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={4}
          >
            <div className="step-number">5</div>
            <div className="step-content">
              <h3>Iterate & Optimize</h3>
              <p>Refine your results through rapid iterations, guided by AI suggestions and your expertise.</p>
            </div>
          </motion.div>
        </div>
      </section>

      {/* Core Services Section */}
      <section className="services-section">
        <motion.h2 
          className="section-title"
          initial="hidden"
          whileInView="visible"
          viewport={{ once: true }}
          variants={fadeInUpVariants}
        >
          Core Services
        </motion.h2>
        
        <div className="services-grid">
          {/* Protein Structure Analysis */}
          <motion.div 
            className="service-card glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={0}
          >
            <div className="service-icon">🧬</div>
            <h3>Protein Structure Analysis</h3>
            <p>Comprehensive tools for analyzing protein structures and predicting their functions.</p>
            
            <div className="capabilities">
              <h4>Key Capabilities:</h4>
              <ul>
                <li>PDB fetching & parsing</li>
                <li>3D visualization</li>
                <li>Secondary structure prediction</li>
                <li>Mutation impact analysis</li>
              </ul>
            </div>
            
            <div className="automation-example">
              <h4>Agent-Powered Automation:</h4>
              <div className="chat-example">
                "Agent, analyze protein X for potential binding pockets."
              </div>
            </div>
            
            <Link to="/protein-analysis" className="service-cta">Explore Protein Analysis Tools</Link>
          </motion.div>
          
          {/* Drug Candidate Identification */}
          <motion.div 
            className="service-card glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={1}
          >
            <div className="service-icon">💊</div>
            <h3>Drug Candidate Identification</h3>
            <p>AI-powered tools to identify and evaluate potential drug candidates efficiently.</p>
            
            <div className="capabilities">
              <h4>Key Capabilities:</h4>
              <ul>
                <li>Virtual screening</li>
                <li>QSAR modeling</li>
                <li>ADMET prediction</li>
                <li>Molecular docking</li>
              </ul>
            </div>
            
            <div className="automation-example">
              <h4>Agent-Powered Automation:</h4>
              <div className="chat-example">
                "Agent, screen compound library against target protein and rank by binding affinity."
              </div>
            </div>
            
            <Link to="/drug-identification" className="service-cta">Explore Drug Identification Tools</Link>
          </motion.div>
          
          {/* Target Optimization */}
          <motion.div 
            className="service-card glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={2}
          >
            <div className="service-icon">🎯</div>
            <h3>Target Optimization</h3>
            <p>Refine and optimize promising drug candidates for improved efficacy and safety.</p>
            
            <div className="capabilities">
              <h4>Key Capabilities:</h4>
              <ul>
                <li>Lead optimization suggestions</li>
                <li>Affinity maturation guidance</li>
                <li>Selectivity profiling</li>
                <li>Resistance mutation analysis</li>
              </ul>
            </div>
            
            <div className="automation-example">
              <h4>Agent-Powered Automation:</h4>
              <div className="chat-example">
                "Agent, suggest R-group modifications to improve compound solubility."
              </div>
            </div>
            
            <Link to="/target-optimization" className="service-cta">Explore Target Optimization Tools</Link>
          </motion.div>
        </div>
      </section>

      {/* Chat Interface Section */}
      <section className="chat-section">
        <motion.div 
          className="chat-content"
          initial="hidden"
          whileInView="visible"
          viewport={{ once: true }}
          variants={fadeInLeftVariants}
        >
          <h2>Converse, Command, Create: Drug Discovery via Intuitive Chat</h2>
          <p>
            Our natural language chat interface serves as your command center, allowing you to:
          </p>
          <ul>
            <li>Initiate complex workflows with simple commands</li>
            <li>Adjust parameters in real-time</li>
            <li>Query results and generate insights</li>
            <li>Launch specialized tools when needed</li>
          </ul>
          <div className="chat-commands">
            <div className="command-example">
              "Compare binding affinities of compounds A and B with target protein."
            </div>
            <div className="command-example">
              "Visualize molecular dynamics simulation and highlight residues 45-60."
            </div>
            <div className="command-example">
              "Generate ADMET predictions for the top 5 candidates."
            </div>
          </div>
        </motion.div>
        <motion.div 
          className="chat-image"
          initial="hidden"
          whileInView="visible"
          viewport={{ once: true }}
          variants={fadeInRightVariants}
        >
          <img src="/images/chatbot.png" alt="Chat Interface" className="chat-interface-img" />
        </motion.div>
      </section>

      {/* Tools & Integrations Section */}
      <section className="tools-section" id="tools">
        <motion.h2 
          className="section-title"
          initial="hidden"
          whileInView="visible"
          viewport={{ once: true }}
          variants={fadeInUpVariants}
        >
          Built on a Foundation of Open Science and Cutting-Edge APIs
        </motion.h2>
        
        <div className="tools-grid">
          <motion.div 
            className="tool-item glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={0}
          >
            <h3>RDKit</h3>
            <p>Powerful cheminformatics and machine learning toolkit</p>
          </motion.div>
          
          <motion.div 
            className="tool-item glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={1}
          >
            <h3>BioPython</h3>
            <p>Tools for biological computation</p>
          </motion.div>
          
          <motion.div 
            className="tool-item glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={2}
          >
            <h3>PDB</h3>
            <p>Protein Data Bank integration</p>
          </motion.div>
          
          <motion.div 
            className="tool-item glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={3}
          >
            <h3>AlphaFold DB</h3>
            <p>Predicted protein structures database</p>
          </motion.div>
          
          <motion.div 
            className="tool-item glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={4}
          >
            <h3>ChEMBL</h3>
            <p>Bioactive molecules database</p>
          </motion.div>
          
          <motion.div 
            className="tool-item glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={5}
          >
            <h3>PubChem</h3>
            <p>Chemical information database</p>
          </motion.div>
          
          <motion.div 
            className="tool-item and-more glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={6}
          >
            <h3>And Many More...</h3>
          </motion.div>
        </div>
        
        <motion.div 
          className="tools-cta"
          initial="hidden"
          whileInView="visible"
          viewport={{ once: true }}
          variants={fadeInUpVariants}
        >
          <Link to="/technology" className="secondary-btn">Learn More About Our Technology Stack</Link>
        </motion.div>
      </section>

      {/* Benefits Section */}
      <section className="benefits-section">
        <motion.h2 
          className="section-title"
          initial="hidden"
          whileInView="visible"
          viewport={{ once: true }}
          variants={fadeInUpVariants}
        >
          Transform Your Drug Discovery Workflow
        </motion.h2>
        
        <div className="benefits-grid">
          <motion.div 
            className="benefit-card glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={0}
          >
            <div className="benefit-icon">⚡</div>
            <h3>Accelerate Research</h3>
            <p>Automate repetitive tasks and streamline workflows to reduce research time by up to 70%.</p>
          </motion.div>
          
          <motion.div 
            className="benefit-card glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={1}
          >
            <div className="benefit-icon">🔍</div>
            <h3>Enhance Precision & Control</h3>
            <p>Maintain complete control through chat-based tweaking and direct access to specialized tools.</p>
          </motion.div>
          
          <motion.div 
            className="benefit-card glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={2}
          >
            <div className="benefit-icon">🧩</div>
            <h3>Simplify Complexity</h3>
            <p>Access sophisticated analysis tools through natural language, eliminating steep learning curves.</p>
          </motion.div>
          
          <motion.div 
            className="benefit-card glassy-feel"
            initial="hidden"
            whileInView="visible"
            viewport={{ once: true }}
            variants={fadeInUpVariants}
            custom={3}
          >
            <div className="benefit-icon">🔄</div>
            <h3>Foster Iteration & Innovation</h3>
            <p>Rapidly test hypotheses and iterate on designs with immediate feedback and AI-powered suggestions.</p>
          </motion.div>
        </div>
      </section>

      {/* Footer */}
      <footer className="main-footer">
        <div className="footer-content">
          <div className="footer-logo">
            <img src="/images/valency_logo_light_600x600.png" alt="Valency Logo" className="footer-logo-img" />
            <span className="footer-logo-text">Valency</span>
          </div>
          
          <div className="footer-links">
            <div className="footer-links-column">
              <h3>Company</h3>
              <Link to="/about">About Us</Link>
              <Link to="/careers">Careers</Link>
              <Link to="/contact">Contact</Link>
            </div>
            
            <div className="footer-links-column">
              <h3>Resources</h3>
              <Link to="/documentation">Documentation</Link>
              <Link to="/tutorials">Tutorials</Link>
              <Link to="/blog">Blog</Link>
            </div>
            
            <div className="footer-links-column">
              <h3>Legal</h3>
              <Link to="/privacy">Privacy Policy</Link>
              <Link to="/terms">Terms of Service</Link>
            </div>
            
            <div className="footer-links-column">
              <h3>Connect</h3>
              <div className="social-links">
                <a href="https://linkedin.com" target="_blank" rel="noopener noreferrer" className="social-link">LinkedIn</a>
                <a href="https://twitter.com" target="_blank" rel="noopener noreferrer" className="social-link">Twitter</a>
                <a href="https://github.com" target="_blank" rel="noopener noreferrer" className="social-link">GitHub</a>
              </div>
            </div>
          </div>
        </div>
        
        <div className="copyright">
          <p>© 2025 Valency. All Rights Reserved.</p>
        </div>
      </footer>
    </div>
  );
}

export default Home;